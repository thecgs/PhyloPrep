import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';

function uleb(bytes, cursor) {
  let value = 0, factor = 1, byte;
  do { byte = bytes[cursor.pos++]; value += (byte & 127) * factor; factor *= 128; } while (byte & 128);
  return value;
}
export function memoryLimits(bytes) {
  const cursor = { pos: 8 };
  while (cursor.pos < bytes.length) {
    const section = bytes[cursor.pos++], size = uleb(bytes, cursor), end = cursor.pos + size;
    if (section === 5) {
      assert.equal(uleb(bytes, cursor), 1, 'Expected one defined memory');
      const flags = uleb(bytes, cursor), initial = uleb(bytes, cursor);
      assert.ok(flags & 1, 'Memory must declare a maximum');
      return { initial, maximum: uleb(bytes, cursor) };
    }
    cursor.pos = end;
  }
  throw new Error('No defined memory section');
}
const manifest = JSON.parse(await readFile('site/aligners/manifest.json', 'utf8'));
assert.equal(manifest.maximumMemoryBytes, 4294967296);
for (const [tool, metadata] of Object.entries(manifest.tools)) {
  for (const file of metadata.files) {
    const bytes = await readFile(`site/aligners/${tool}/${metadata.packageVersion}/${file.name}`);
    assert.equal(createHash('sha256').update(bytes).digest('hex'), file.sha256);
    if (file.name.endsWith('.wasm')) {
      assert.equal(memoryLimits(bytes).maximum, 65536, `${tool}/${file.name} must allow 4 GiB`);
      assert.ok(WebAssembly.validate(bytes), `Invalid WebAssembly: ${file.name}`);
    }
  }
}
console.log('All three compiled aligners declare a 4 GiB maximum and match their checksums.');
