const CACHE_NAME = "phyloprep-web-v76";
const APP_SHELL = [
  "./",
  "./index.html",
  "./styles.css",
  "./app.js",
  "./browser-inputs.js",
  "./browser-downloads.js",
  "./browser-results.js",
  "./browser-save.js",
  "./browser-hash.js",
  "./browser-reuse.js",
  "./orthogroup-defaults.js",
  "./pyodide-worker.js",
  "./browser-disk-fs.js",
  "./alignment-worker.js",
  "./manifest.webmanifest",
  "./icon.svg",
];
const PYODIDE_PREFIX = "https://cdn.jsdelivr.net/pyodide/v0.27.2/full/";

self.addEventListener("install", (event) => {
  event.waitUntil(caches.open(CACHE_NAME).then((cache) => cache.addAll(APP_SHELL)).then(() => self.skipWaiting()));
});

self.addEventListener("activate", (event) => {
  event.waitUntil(caches.keys().then((keys) => Promise.all(keys.filter((key) => key.startsWith("phyloprep-web-") && key !== CACHE_NAME).map((key) => caches.delete(key)))).then(() => self.clients.claim()));
});

self.addEventListener("fetch", (event) => {
  if (event.request.method !== "GET") return;
  const url = new URL(event.request.url);
  const cacheable = url.origin === self.location.origin || url.href.startsWith(PYODIDE_PREFIX);
  if (!cacheable) return;
  if (event.request.cache === "no-store") {
    event.respondWith(fetch(event.request));
    return;
  }
  const refresh = event.request.cache === "no-cache" || url.pathname.endsWith("/engine/manifest.json");
  const network = () => fetch(event.request).then(async response => {
    if (response.ok || response.type === "opaque") {
      const cache = await caches.open(CACHE_NAME);
      await cache.put(event.request, response.clone());
    }
    return response;
  });
  event.respondWith((async () => {
    const cached = await caches.match(event.request);
    if (refresh) {
      try {
        const response = await network();
        if (response.ok || !cached) return response;
      } catch (error) { if (!cached) throw error; }
      return cached;
    }
    if (cached) return cached;
    return network();
  })());
});
