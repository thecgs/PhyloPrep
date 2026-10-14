// Native browser download helpers. No timeout can reliably identify completion.
(() => {
  const urls = new WeakMap();
  function requestDownload(makeBlob, name, owner) {
    let url = urls.get(owner);
    if (!url) {
      url = URL.createObjectURL(makeBlob());
      urls.set(owner, url);
    }
    const link = document.createElement('a');
    link.href = url;
    link.download = name;
    link.hidden = true;
    document.body.append(link);
    try { link.click(); }
    finally { link.remove(); }
    // Keep the URL valid for retries and delayed save dialogs. The browser
    // releases it when this document is unloaded; repeated clicks reuse it.
  }
  let table;
  function crc32(bytes) {
    if (!table) table = Array.from({ length: 256 }, (_, index) => {
      let value = index;
      for (let bit = 0; bit < 8; bit++) value = value & 1 ? (value >>> 1) ^ 0xedb88320 : value >>> 1;
      return value >>> 0;
    });
    let value = 0xffffffff;
    for (const byte of bytes) value = (value >>> 8) ^ table[(value ^ byte) & 0xff];
    return (value ^ 0xffffffff) >>> 0;
  }
  function checkZipCount(files) {
    if (files.length >= 0xffff) throw new Error('This archive exceeds the ZIP32 file-count limit. Select an output folder and save all files there.');
  }
  function* zipParts(files) {
    checkZipCount(files);
    const encoder = new TextEncoder(), directory = [];
    let offset = 0;
    for (const file of files) {
      const name = encoder.encode(file.name), content = encoder.encode(file.content);
      if (name.length > 0xffff || content.length >= 0xffffffff || offset + 30 + name.length + content.length >= 0xffffffff) {
        throw new Error('This archive exceeds ZIP32 limits. Select an output folder and save all files there.');
      }
      const checksum = crc32(content);
      const local = new Uint8Array(30 + name.length), view = new DataView(local.buffer);
      view.setUint32(0, 0x04034b50, true); view.setUint16(4, 20, true);
      view.setUint16(6, 0x0800, true); // UTF-8 names
      view.setUint32(14, checksum, true); view.setUint32(18, content.length, true);
      view.setUint32(22, content.length, true); view.setUint16(26, name.length, true);
      local.set(name, 30);
      yield local; yield content;
      const central = new Uint8Array(46 + name.length), cv = new DataView(central.buffer);
      cv.setUint32(0, 0x02014b50, true); cv.setUint16(4, 20, true); cv.setUint16(6, 20, true);
      cv.setUint16(8, 0x0800, true); cv.setUint32(16, checksum, true);
      cv.setUint32(20, content.length, true); cv.setUint32(24, content.length, true);
      cv.setUint16(28, name.length, true); cv.setUint32(42, offset, true);
      central.set(name, 46); directory.push(central);
      offset += local.length + content.length;
    }
    const centralSize = directory.reduce((sum, part) => sum + part.length, 0);
    if (offset + centralSize + 22 >= 0xffffffff) throw new Error('This archive exceeds ZIP32 size limits. Save files to an output folder.');
    for (const part of directory) yield part;
    const end = new Uint8Array(22), ev = new DataView(end.buffer);
    ev.setUint32(0, 0x06054b50, true); ev.setUint16(8, files.length, true); ev.setUint16(10, files.length, true);
    ev.setUint32(12, centralSize, true); ev.setUint32(16, offset, true);
    yield end;
  }
  async function writeZip(files, writable) {
    try {
      for (const part of zipParts(files)) await writable.write(part);
      await writable.close();
    } catch (error) {
      try { await writable.abort(); } catch {}
      throw error;
    }
  }
  globalThis.PhyloPrepDownloads = { requestDownload, zipParts, writeZip, checkZipCount };
})();
