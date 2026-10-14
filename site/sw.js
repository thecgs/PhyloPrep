const CACHE_NAME = "phyloprep-web-v30";
const APP_SHELL = [
  "./",
  "./index.html",
  "./styles.css",
  "./app.js",
  "./pyodide-worker.js",
  "./alignment-worker.js",
  "./browser_runner.py",
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
  event.respondWith(caches.match(event.request).then((cached) => {
    const network = fetch(event.request).then((response) => {
      if (response.ok || response.type === "opaque") caches.open(CACHE_NAME).then((cache) => cache.put(event.request, response.clone()));
      return response;
    });
    return cached || network;
  }));
});
