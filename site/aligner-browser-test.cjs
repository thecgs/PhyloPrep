// Executed after building the aligners. Playwright is a CI-only test dependency.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { chromium } = require('playwright');
const root = path.resolve('site');
const server = http.createServer((req, res) => {
  const file = path.resolve(root, '.' + new URL(req.url, 'http://localhost').pathname);
  if (!file.startsWith(root + path.sep)) { res.writeHead(403).end(); return; }
  if (!fs.existsSync(file) || !fs.statSync(file).isFile()) { res.writeHead(404).end(); return; }
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Content-Type', file.endsWith('.wasm') ? 'application/wasm' : file.endsWith('.js') ? 'application/javascript' : 'text/html');
  fs.createReadStream(file).pipe(res);
});

(async () => {
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const origin = `http://127.0.0.1:${server.address().port}`;
  const proxyUrl = process.env.PHYLOPREP_TEST_PROXY ? new URL(process.env.PHYLOPREP_TEST_PROXY) : null;
  const proxy = proxyUrl ? { server: proxyUrl.origin, username: decodeURIComponent(proxyUrl.username),
    password: decodeURIComponent(proxyUrl.password), bypass: '127.0.0.1,localhost' } : undefined;
  const browser = await chromium.launch({ headless: true, args: ['--no-sandbox'], proxy });
  try {
    const context = await browser.newContext();
    context.on('requestfailed', request => console.error('Request failed:', request.url(), request.failure()?.errorText));
    // Optional byte-for-byte CDN mirror for environments whose browser cannot
    // use the host's network proxy. Production workers and wasm stay unchanged.
    if (process.env.PHYLOPREP_TEST_VENDOR_DIR) {
      await context.route('https://biowasm.com/**', async route => {
        const basename = path.basename(new URL(route.request().url()).pathname);
        if (!['aioli.js', 'base.js', 'base.wasm'].includes(basename)) throw new Error(`Unexpected CDN asset ${basename}`);
        await route.fulfill({path: path.join(process.env.PHYLOPREP_TEST_VENDOR_DIR, basename),
          contentType: basename.endsWith('.wasm') ? 'application/wasm' : 'application/javascript'});
      });
    }
    const page = await context.newPage();
    page.on('console', message => { if (message.type() === 'error') console.error(message.text()); });
    await page.goto(`${origin}/index.html`);
    // The memory probe uses the same compiler and ceiling and touches >2 GiB.
    const probe = await page.evaluate(async origin => {
      const source = `importScripts('https://biowasm.com/cdn/v3/aioli.js');
        (async () => { try {
          const cli = await new Aioli([{tool:'probe',version:'1',program:'probe',urlPrefix:'${origin}/aligners/probe/1'}], {printInterleaved:false});
          const result = await cli.exec('probe'); await cli.close(); postMessage(result);
        } catch(e) { postMessage({error:String(e)}); } })();`;
      return await new Promise((resolve, reject) => {
        const url = URL.createObjectURL(new Blob([source], {type:'application/javascript'}));
        const worker = new Worker(url);
        const timer = setTimeout(() => { worker.terminate(); reject(new Error('Memory probe timed out')); }, 180000);
        worker.onmessage = e => { clearTimeout(timer); worker.terminate(); URL.revokeObjectURL(url); resolve(e.data); };
        worker.onerror = e => { clearTimeout(timer); worker.terminate(); reject(new Error(e.message)); };
      });
    }, origin);
    assert.ok(!probe.error, probe.error);
    assert.match(probe.stdout, /END=29/);
    const heap = Number(/HEAP_BYTES=(\d+)/.exec(probe.stdout)?.[1]);
    assert.ok(heap > 2147483648 && heap <= 4294967296, `Unexpected probe heap: ${heap}`);
    console.log(`Browser memory probe crossed 2 GiB: ${heap} bytes`);
    // Run real production workers twice to check tool loading, output, and reuse.
    const input = '>a\nACDEFGHIKLMNPQRSTVWY\n>b\nACDEFGHIKLMNPQRSAVWY\n>c\nACDEFGHIKLMNPQRSTVW\n';
    for (const tool of ['mafft', 'muscle']) {
      const outputs = await page.evaluate(async ({tool, input}) => {
        const worker = new Worker('alignment-worker.js');
        try {
          const run = runId => new Promise((resolve, reject) => {
            const timer = setTimeout(() => reject(new Error('Alignment timed out')), 180000);
            worker.onmessage = e => { clearTimeout(timer); e.data.type === 'error' ? reject(new Error(e.data.message)) : resolve(e.data.alignment); };
            worker.onerror = e => { clearTimeout(timer); reject(new Error(e.message)); };
            worker.postMessage({type:'align',runId,tool,seqtype:'prot',sequenceText:input});
          });
          return [await run(1), await run(2)];
        } finally { worker.terminate(); }
      }, {tool, input});
      for (const output of outputs) {
        const records = output.trim().split(/^>/m).slice(1).map(s => s.split(/\r?\n/).slice(1).join('').replace(/\s/g,''));
        assert.equal(records.length, 3, `${tool} lost taxa`);
        assert.equal(new Set(records.map(s => s.length)).size, 1, `${tool} produced unequal lengths`);
        assert.deepEqual(records.map(s => s.replace(/-/g, '')).sort(), ['ACDEFGHIKLMNPQRSTVWY','ACDEFGHIKLMNPQRSAVWY','ACDEFGHIKLMNPQRSTVW'].sort());
      }
      console.log(`${tool}: two real browser alignments passed`);
    }
  } finally { await browser.close(); server.close(); }
})().catch(error => { console.error(error); server.close(); process.exitCode = 1; });
