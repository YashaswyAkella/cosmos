/* ===========================================================
   Our Universe — 3D engine + interface
   A small perspective renderer on canvas 2D: real Keplerian
   orbits, true 3D star positions, an orbiting camera that
   spans from the surface of the Sun to the galactic centre.
   =========================================================== */
(function () {
'use strict';

/* ============================================================
   1. Math & coordinates
   ============================================================ */
const DEG = Math.PI / 180;
const TAU = Math.PI * 2;
const OBLIQ = 23.4392911 * DEG;          // obliquity of the ecliptic, J2000
const R_SUN_AU = 0.00465047;             // solar radius in AU
const J2000 = 2451545.0;

const clamp = (v, a, b) => v < a ? a : (v > b ? b : v);
const lerp = (a, b, t) => a + (b - a) * t;

function jdFromDate(d) { return d.getTime() / 86400000 + 2440587.5; }
function dateFromJd(jd) { return new Date((jd - 2440587.5) * 86400000); }

/* equatorial cartesian -> ecliptic cartesian */
function eqToEcl(x, y, z) {
  const c = Math.cos(OBLIQ), s = Math.sin(OBLIQ);
  return { x: x, y: y * c + z * s, z: -y * s + z * c };
}
/* RA (hours) / Dec (deg) / distance -> ecliptic cartesian */
function raDecToVec(raH, decDeg, dist) {
  const ra = raH * 15 * DEG, dec = decDeg * DEG;
  const cd = Math.cos(dec);
  return eqToEcl(dist * cd * Math.cos(ra), dist * cd * Math.sin(ra), dist * Math.sin(dec));
}

/* galactic -> equatorial (Hipparcos rotation, transposed) then -> ecliptic */
const GAL = [
  [-0.0548755604, +0.4941094279, -0.8676661490],
  [-0.8734370902, -0.4448296300, -0.1980763734],
  [-0.4838350155, +0.7469822445, +0.4559837762]
];
function galacticToEcl(l, b) {
  const cb = Math.cos(b);
  const gx = cb * Math.cos(l), gy = cb * Math.sin(l), gz = Math.sin(b);
  const ex = GAL[0][0] * gx + GAL[0][1] * gy + GAL[0][2] * gz;
  const ey = GAL[1][0] * gx + GAL[1][1] * gy + GAL[1][2] * gz;
  const ez = GAL[2][0] * gx + GAL[2][1] * gy + GAL[2][2] * gz;
  return eqToEcl(ex, ey, ez);
}

function norm(v) {
  const m = Math.hypot(v.x, v.y, v.z) || 1;
  return { x: v.x / m, y: v.y / m, z: v.z / m };
}
function cross(a, b) {
  return { x: a.y * b.z - a.z * b.y, y: a.z * b.x - a.x * b.z, z: a.x * b.y - a.y * b.x };
}

/* Kepler: elements -> heliocentric ecliptic position (AU) */
function keplerPos(el, T) {
  const a = el.a + (el.aD || 0) * T;
  const e = el.e + (el.eD || 0) * T;
  const I = (el.i + (el.iD || 0) * T) * DEG;
  const L = (el.L + (el.LD || 0) * T) * DEG;
  const w = (el.w + (el.wD || 0) * T) * DEG;
  const O = (el.O + (el.OD || 0) * T) * DEG;

  let M = L - w;
  M = M - TAU * Math.floor((M + Math.PI) / TAU);

  let E = M + e * Math.sin(M);
  for (let k = 0; k < 12; k++) {
    const d = (M - (E - e * Math.sin(E))) / (1 - e * Math.cos(E));
    E += d;
    if (Math.abs(d) < 1e-12) break;
  }
  return orbitalToEcl(a, e, E, w - O, O, I);
}
/* position on the ellipse for a given eccentric anomaly */
function orbitalToEcl(a, e, E, argPeri, O, I) {
  const xp = a * (Math.cos(E) - e);
  const yp = a * Math.sqrt(1 - e * e) * Math.sin(E);
  const cw = Math.cos(argPeri), sw = Math.sin(argPeri);
  const cO = Math.cos(O), sO = Math.sin(O);
  const cI = Math.cos(I), sI = Math.sin(I);
  return {
    x: (cw * cO - sw * sO * cI) * xp + (-sw * cO - cw * sO * cI) * yp,
    y: (cw * sO + sw * cO * cI) * xp + (-sw * sO + cw * cO * cI) * yp,
    z: (sw * sI) * xp + (cw * sI) * yp
  };
}

function mulberry32(seed) {
  return function () {
    seed |= 0; seed = seed + 0x6D2B79F5 | 0;
    let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

/* ============================================================
   2. Formatting
   ============================================================ */
function fmtDist(au) {
  if (!isFinite(au)) return '—';
  if (au < 0.0007) return Math.round(au * AU_KM).toLocaleString() + ' km';
  if (au < 0.02) return (au * AU_KM / 1000).toFixed(1).replace(/\.0$/, '') + ' thousand km';
  if (au < 1200) return au.toFixed(au < 10 ? 3 : 2) + ' AU';
  const ly = au / LY_AU;
  if (ly < 1000) return ly.toFixed(ly < 10 ? 3 : 1) + ' light years';
  if (ly < 1e6) return Math.round(ly).toLocaleString() + ' light years';
  if (ly < 1e9) return (ly / 1e6).toFixed(ly < 1e7 ? 2 : ly < 1e8 ? 1 : 0) + ' million light years';
  return (ly / 1e9).toFixed(1) + ' billion light years';
}
function fmtDate(jd) {
  const d = dateFromJd(jd);
  return d.toLocaleDateString(undefined, { day: '2-digit', month: 'short', year: 'numeric', timeZone: 'UTC' }) +
         ' · ' + String(d.getUTCHours()).padStart(2, '0') + ':' + String(d.getUTCMinutes()).padStart(2, '0');
}

/* ============================================================
   3. Build the scene graph from the catalog
   ============================================================ */
const SPEC = {
  O: { c: '#9bb0ff', r: 8.0 }, B: { c: '#aabfff', r: 4.0 }, A: { c: '#cad7ff', r: 1.8 },
  F: { c: '#f8f7ff', r: 1.3 }, G: { c: '#fff4ea', r: 1.0 }, K: { c: '#ffd2a1', r: 0.78 },
  M: { c: '#ff9d6c', r: 0.35 }, D: { c: '#dfe8ff', r: 0.013 }
};

const objects = [];
const byId = {};

function addObject(o) {
  o.pos = { x: 0, y: 0, z: 0 };
  objects.push(o);
  byId[o.id] = o;
  return o;
}
function poleVec(raDeg, decDeg) {
  if (raDeg == null) return { x: 0, y: 0, z: 1 };
  return norm(raDecToVec(raDeg / 15, decDeg, 1));
}

/* --- Sun --- */
addObject({
  id: SUN.id, name: SUN.name, type: SUN.type, kind: 'sun', category: 'solarsystem',
  color: SUN.color, glow: SUN.glow, radiusAU: SUN.radiusKm / AU_KM,
  pole: poleVec(SUN.poleRA, SUN.poleDec), data: SUN, prio: 100, minPx: 3.5
});

/* --- Planets, dwarf planets and their moons --- */
function addPlanetLike(p, category) {
  const o = addObject({
    id: p.id, name: p.name, type: p.type, kind: p.kind, category: category,
    color: p.color, color2: p.color2 || p.color, radiusAU: p.radiusKm / AU_KM,
    el: p.el, pole: poleVec(p.poleRA, p.poleDec), rings: p.rings, bands: p.bands,
    data: p, prio: p.kind === 'planet' ? 90 : 72, minPx: p.kind === 'planet' ? 2.6 : 2.1
  });
  (p.moons || []).forEach(function (m) {
    addObject({
      id: m.id, name: m.name, type: 'Moon of ' + p.name, kind: 'moon', category: 'moons',
      color: m.color, color2: m.color, radiusAU: m.radiusKm / AU_KM, pole: m.poleRA != null ? poleVec(m.poleRA, m.poleDec) : null, tidal: !!m.tidal,
      parent: p.id, aAU: m.aKm / AU_KM, periodDays: m.periodDays, inc: (m.inc || 0) * DEG,
      refPlane: m.refPlane,
      data: m, prio: 48, minPx: 1.8
    });
  });
  return o;
}
PLANETS.forEach(function (p) { addPlanetLike(p, 'solarsystem'); });
DWARFS.forEach(function (p) { addPlanetLike(p, 'dwarfs'); });

/* --- Spacecraft --- */
PROBES.forEach(function (p) {
  const o = {
    id: p.id, name: p.name, type: p.type, kind: 'probe', category: p.cat || 'probes',
    color: p.color, radiusAU: 0, data: p, prio: 64, minPx: 0, mode: p.mode
  };
  if (p.mode === 'escape') {
    o.dir = norm(raDecToVec(p.ra, p.dec, 1));
    o.epochJd = jdFromDate(new Date(p.epoch + 'T00:00:00Z'));
    o.distAU = p.distAU;
    o.speedAUday = p.speedAUyr / 365.25;
  } else if (p.mode === 'kepler') {
    o.el = p.el;
  } else if (p.mode === 'earthorbit' || p.mode === 'orbiter') {
    o.aAU = p.aKm / AU_KM;
    o.periodDays = p.periodDays;
    o.inc = (p.inc || 0) * DEG;
    o.ecc = p.e || 0;
    o.argPeri = (p.argPeri || 0) * DEG;
    o.parent = p.parent || 'earth';
  } else if (p.mode === 'atbody') {
    o.parent = p.parent;
  } else if (p.mode === 'trailing') {
    o.launchJd = jdFromDate(new Date(p.launch + 'T00:00:00Z'));
    o.driftDegYr = p.driftDegYr;
  } else if (p.mode === 'l2' || p.mode === 'l2halo') {
    o.offsetAU = p.offsetAU;
    o.haloAU = p.haloAU;
    o.haloDays = p.haloDays;
    if (p.launch) {
      o.launchJd = jdFromDate(new Date(p.launch + 'T00:00:00Z'));
      o.arriveJd = jdFromDate(new Date(p.arrive + 'T00:00:00Z'));
    }
  }
  addObject(o);
});

/* --- Stars --- */
STARS.forEach(function (s) {
  const sp = SPEC[s.spec] || SPEC.G;
  const distAU = s.d * LY_AU;
  const v = raDecToVec(s.ra, s.dec, distAU);
  const o = addObject({
    id: s.id, name: s.name, type: (s.spec === 'D' ? 'White dwarf' : 'Star') + (s.note ? ' · ' + s.note : ''),
    kind: 'star', category: 'stars', color: sp.c,
    radiusAU: (s.rSun != null ? s.rSun : sp.r) * R_SUN_AU,
    absMag: s.mag - 5 * Math.log10(s.d / PC_LY) + 5, appMag: s.mag,
    fixed: v, data: s, prio: 40 - s.mag, minPx: 1.2
  });
  o.pos = { x: v.x, y: v.y, z: v.z };
});

/* --- Black holes --- */
EXOTIC.forEach(function (b) {
  const distAU = b.d * LY_AU;
  const v = raDecToVec(b.ra, b.dec, distAU);
  const o = addObject({
    id: b.id, name: b.name, type: b.type, kind: 'blackhole', category: 'blackholes',
    color: b.color, radiusAU: (b.shadowKm / 2) / AU_KM, fixed: v, data: b, prio: 95, minPx: 2
  });
  o.pos = { x: v.x, y: v.y, z: v.z };
});

/* --- Stellar remnants: what stars leave behind --- */
if (typeof REMNANTS !== 'undefined') REMNANTS.forEach(function (r) {
  const v = raDecToVec(r.ra, r.dec, r.d * LY_AU);

  if (r.sub === 'stellarbh') {
    /* Schwarzschild radius is 2GM/c^2 = 2.953 km per solar mass, so the horizon is
       5.906 km across per solar mass — a few dozen km, not millions like Sgr A* */
    const horizonKm = 5.906 * r.mass;
    const o = addObject({
      id: r.id, name: r.name, type: 'Stellar-mass black hole', kind: 'blackhole',
      category: 'blackholes', color: r.color, radiusAU: (horizonKm / 2) / AU_KM,
      fixed: v, data: r, prio: 66, minPx: 2
    });
    o.pos = o.fixed;
    r.stats = Object.assign({ 'Distance from Sun': fmtLy(r.d) }, r.stats || {});
    return;
  }

  if (r.sub === 'whitedwarf') {
    const am = r.absMag != null ? r.absMag : 12;
    const o = addObject({
      id: r.id, name: r.name, type: 'White dwarf', kind: 'star', category: 'remnants',
      color: r.color, radiusAU: r.rKm / AU_KM,
      absMag: am, appMag: am - 5 * Math.log10(r.d / PC_LY) + 5,
      fixed: v, data: r, prio: 44, minPx: 1.4
    });
    o.pos = o.fixed;
    r.stats = Object.assign({ 'Distance from Sun': fmtLy(r.d) }, r.stats || {});
    return;
  }

  /* neutron stars, pulsars and magnetars: about 22 km across, self-luminous */
  const TYPE = { pulsar: 'Pulsar', neutron: 'Neutron star', magnetar: 'Magnetar' };
  const o = addObject({
    id: r.id, name: r.name, type: TYPE[r.sub] || 'Neutron star', kind: 'neutron',
    sub: r.sub, category: 'remnants', color: r.color,
    radiusAU: r.rKm / AU_KM, beams: !!r.beams,
    fixed: v, data: r, prio: 68, minPx: 1.8
  });
  o.pos = o.fixed;
  /* deterministic spin axis plus a perpendicular basis, for the beams */
  let seed = 0;
  for (let i = 0; i < r.id.length; i++) seed += r.id.charCodeAt(i) * (i + 3);
  const rng = mulberry32(seed);
  const cz = rng() * 1.7 - 0.85, ph = rng() * TAU, sp = Math.sqrt(1 - cz * cz);
  o.spinAxis = norm({ x: sp * Math.cos(ph), y: sp * Math.sin(ph), z: cz });
  let u = cross({ x: 0, y: 0, z: 1 }, o.spinAxis);
  if (Math.hypot(u.x, u.y, u.z) < 1e-6) u = { x: 1, y: 0, z: 0 };
  o.spinU = norm(u);
  o.spinV = cross(o.spinAxis, o.spinU);

  const st = { 'Distance from Sun': fmtLy(r.d) };
  if (r.spinS) st['Rotation period'] = r.spinS < 1 ? (r.spinS * 1000).toFixed(2) + ' ms'
                                                  : r.spinS.toFixed(2) + ' s';
  if (r.rKm) st['Radius'] = r.rKm + ' km';
  r.stats = Object.assign(st, r.stats || {});
});

/* --- Deep sky: nebulae, clusters, galaxies --- */
const DS_TYPE = {
  spiral: 'Spiral galaxy', galaxy: 'Galaxy', nebula: 'Nebula',
  planetary: 'Planetary nebula', supernova: 'Supernova remnant',
  globular: 'Globular cluster', open: 'Open cluster'
};
function fmtLy(ly) {
  if (ly >= 1e6) return (ly / 1e6).toFixed(2).replace(/\.?0+$/, '') + ' million ly';
  return Math.round(ly).toLocaleString() + ' ly';
}
/* Most deep-sky facts follow from the catalogue row, so build the stat
   table rather than repeating it 110 times in the data file. */
function deepSkyStats(d) {
  const st = { 'Type': DS_TYPE[d.sub] || 'Deep sky object' };
  if (d.con) st['Constellation'] = d.con;
  st['Distance'] = fmtLy(d.d);
  if (d.arcmin) st['Apparent size'] = d.arcmin >= 60 ? (d.arcmin / 60).toFixed(1) + '°' : d.arcmin + '′';
  if (d.mag != null) st['Magnitude'] = d.mag.toFixed(1);
  const extra = d.stats || {};
  Object.keys(extra).forEach(function (k) { st[k] = extra[k]; });
  return st;
}

DEEP_SKY.forEach(function (d) {
  const v = raDecToVec(d.ra, d.dec, d.d * LY_AU);
  /* physical size from angular size where that is what the catalogue gives */
  const sizeLy = d.size != null ? d.size : d.d * (d.arcmin / 60) * DEG;
  d.stats = deepSkyStats(d);
  const o = addObject({
    id: d.id, name: d.name, type: DS_TYPE[d.sub] || 'Deep sky object',
    kind: 'deepsky', sub: d.sub, category: d.cat,
    color: d.color, color2: d.color2 || d.color,
    radiusAU: (sizeLy / 2) * LY_AU,
    fixed: v, data: d, prio: 56, minPx: 2.2
  });
  o.pos = { x: v.x, y: v.y, z: v.z };

  const rng = mulberry32(d.id.length * 7919 + d.name.length * 131 + 7);
  if (d.sub === 'globular' || d.sub === 'open') {
    /* member stars, in unit coordinates scaled by the screen radius */
    const n = d.sub === 'globular' ? 110 : 30;
    o.members = [];
    for (let k = 0; k < n; k++) {
      const rr = d.sub === 'globular' ? Math.pow(rng(), 2.3) : Math.pow(rng(), 1.5);
      const cz = rng() * 2 - 1, ph = rng() * TAU, sp = Math.sqrt(1 - cz * cz);
      o.members.push({ x: sp * Math.cos(ph) * rr, y: sp * Math.sin(ph) * rr, b: 0.3 + rng() * 0.7 });
    }
  } else if (d.sub !== 'planetary') {
    /* a nebula is a handful of overlapping puffs */
    o.puffs = [];
    for (let k = 0; k < 7; k++) {
      o.puffs.push({
        x: (rng() - 0.5) * 1.05, y: (rng() - 0.5) * 1.05,
        r: 0.30 + rng() * 0.45, c: k % 2 ? (d.color2 || d.color) : d.color,
        a: 0.20 + rng() * 0.22
      });
    }
  }
});

/* --- Exoplanets: real systems, orbiting their real host stars --- */
if (typeof EXOPLANETS !== 'undefined') EXOPLANETS.forEach(function (e) {
  const host = byId[e.host];
  if (!host) return;
  /* one shared orbital plane per system, deterministic from the host id */
  let seed = 0;
  for (let i = 0; i < e.host.length; i++) seed += e.host.charCodeAt(i) * (i + 7);
  const rng = mulberry32(seed);
  const cz = rng() * 1.5 - 0.75, ph = rng() * TAU, sp = Math.sqrt(1 - cz * cz);
  const pole = norm({ x: sp * Math.cos(ph), y: sp * Math.sin(ph), z: cz });
  let u = cross({ x: 0, y: 0, z: 1 }, pole);
  if (Math.hypot(u.x, u.y, u.z) < 1e-6) u = { x: 1, y: 0, z: 0 };
  u = norm(u);
  const v = cross(pole, u);

  const o = addObject({
    id: e.id, name: e.name, type: 'Exoplanet · ' + host.name, kind: 'exoplanet',
    category: 'exoplanets', color: e.color, color2: e.color,
    radiusAU: (e.rE * 6371) / AU_KM,
    parent: e.host, lightId: e.host, aAU: e.a, periodDays: e.period,
    ou: u, ov: v, data: e, prio: 62, minPx: 2.0
  });
  e.stats = Object.assign({
    'Host star': host.name,
    'Distance': fmtLy(host.data.d),
    'Orbit radius': e.a < 0.1 ? (e.a).toFixed(4) + ' AU' : (e.a).toFixed(3) + ' AU',
    'Radius': e.rE.toFixed(2) + ' × Earth'
  }, e.stats || {});
});

/* --- Deep fields: tiny patches of sky that were stared at for days --- */
const deepFieldList = [];
if (typeof DEEP_FIELDS !== 'undefined') DEEP_FIELDS.forEach(function (f) {
  const dir = norm(raDecToVec(f.ra, f.dec, 1));
  f.scopes = f.scopes || [f.scope];
  f.scope = f.scopes[0];
  const wide = f.arcmin >= 90 ? (f.arcmin / 60).toFixed(1) + '°' : f.arcmin + ' arcmin';
  const size = f.arcminH ? f.arcmin + ' × ' + f.arcminH + ' arcmin' : wide;
  f.stats = Object.assign({ 'Constellation': f.con, 'Seen by': f.scopes.join(', ').replace('JWST', 'Webb'), 'Field width': size }, f.stats || {});
  (f.highlights || []).forEach(function (h) { if (h.ra != null) h.dir = norm(raDecToVec(h.ra, h.dec, 1)); });
  const label = f.kind === 'target' ? 'Telescope view · ' + f.year : f.planned ? 'Planned survey · ' + f.year : 'Deep field · ' + f.year;
  const o = addObject({
    id: f.id, name: f.name, type: label, kind: 'deepfield',
    category: 'deepfields', color: f.color, radiusAU: 0, dir: dir, angArcmin: Math.max(f.arcmin, f.arcminH || 0),
    fixed: { x: dir.x * 1e11, y: dir.y * 1e11, z: dir.z * 1e11 },
    data: f, prio: 30, minPx: 0, noRender: true
  });
  o.pos = o.fixed;
  deepFieldList.push(o);
});

/* --- Constellations: names only, no stick figures --- */
const constellationList = [];
if (typeof CONSTELLATIONS !== 'undefined') {
  const ranked = CONSTELLATIONS.slice().sort(function (a, b) { return b.area - a.area; });
  CONSTELLATIONS.forEach(function (c) {
    const dir = norm(raDecToVec(c.ra, c.dec, 1));
    const rank = ranked.indexOf(c) + 1;
    c.stats = {
      'Area': c.area.toLocaleString() + ' sq°',
      'Size rank': rank + ' of 88',
      'Sky': c.dec >= 0 ? 'Northern' : 'Southern',
      'Centre': c.ra.toFixed(1) + 'h, ' + c.dec.toFixed(1) + '°'
    };
    if (!c.desc) {
      c.desc = 'One of the 88 constellations, covering ' + c.area.toLocaleString() +
               ' square degrees of the ' + (c.dec >= 0 ? 'northern' : 'southern') + ' sky.';
    }
    const o = addObject({
      id: c.id, name: c.name, type: 'Constellation', kind: 'constellation',
      category: 'constellations', color: '#95a6d4', radiusAU: 0,
      dir: dir, label: c.name.toUpperCase(),
      /* bigger figures get placed first, with a nudge for the famous ones */
      weight: c.area + (c.key ? 900 : 0),
      fixed: { x: dir.x * 1e11, y: dir.y * 1e11, z: dir.z * 1e11 },
      data: c, prio: 20, minPx: 0, noRender: true
    });
    o.pos = o.fixed;
    constellationList.push(o);
  });
  constellationList.sort(function (a, b) { return b.weight - a.weight; });
}

/* --- Belts (particle fields) --- */
const beltFields = [];
Object.keys(BELTS).forEach(function (key) {
  const cfg = BELTS[key];
  const rng = mulberry32(key.length * 9127 + 17);
  const n = cfg.count;
  const f = {
    cfg: cfg, n: n,
    r: new Float64Array(n), th: new Float64Array(n), w: new Float64Array(n),
    ux: new Float64Array(n), uy: new Float64Array(n), uz: new Float64Array(n),
    vx: new Float64Array(n), vy: new Float64Array(n), vz: new Float64Array(n),
    b: new Float64Array(n)
  };
  for (let k = 0; k < n; k++) {
    const r = cfg.rMin + (cfg.rMax - cfg.rMin) * (cfg.spherical ? Math.pow(rng(), 0.6) : rng());
    f.r[k] = r;
    f.th[k] = rng() * TAU;
    f.b[k] = 0.25 + rng() * 0.75;
    if (cfg.spherical) {
      /* isotropic shell: a fixed random direction, essentially static */
      const u = rng() * 2 - 1, ph = rng() * TAU, s = Math.sqrt(1 - u * u);
      f.ux[k] = s * Math.cos(ph); f.uy[k] = s * Math.sin(ph); f.uz[k] = u;
      f.vx[k] = 0; f.vy[k] = 0; f.vz[k] = 0; f.w[k] = 0;
    } else {
      const node = rng() * TAU;
      const inc = (rng() + rng() + rng() - 1.5) * cfg.incSigma * DEG;
      const cO = Math.cos(node), sO = Math.sin(node), ci = Math.cos(inc), si = Math.sin(inc);
      f.ux[k] = cO; f.uy[k] = sO; f.uz[k] = 0;
      f.vx[k] = -sO * ci; f.vy[k] = cO * ci; f.vz[k] = si;
      f.w[k] = TAU / (365.25 * Math.pow(r, 1.5));    // Kepler's third law
    }
  }
  beltFields.push(f);
  /* a selectable entry so belts show up in search and the browse panel */
  addObject({
    id: cfg.id, name: cfg.name, type: 'Region', kind: 'region', category: 'belts',
    color: cfg.color, radiusAU: 0, regionRadius: (cfg.rMin + cfg.rMax) / 2,
    data: { desc: cfg.desc, stats: cfg.stats }, prio: 30, minPx: 0, noRender: true
  });
});


/* --- The Milky Way as a whole: a selectable object, and a model of the galaxy
   for when the camera is far enough out to see it from outside. The model is
   an impression: arm and bar geometry follow survey measurements (Reid et al.
   2019 for the arms, a 3.5 kpc bar at 30° to the Sun-centre line), the light is
   procedural. Rotation is clockwise seen from the north galactic pole, so the
   trailing arms wind outward counter-clockwise in our frame. --- */
const MW = (function () {
  if (typeof MILKY_WAY === 'undefined' || !byId.sgrA) return null;
  const c = byId.sgrA.fixed || byId.sgrA.pos;            /* the centre is the black hole */
  const KPC = 1000 * PC_LY * LY_AU;                         /* one kiloparsec in AU */
  const ex = norm({ x: -c.x, y: -c.y, z: -c.z });           /* centre -> Sun */
  let ez = galacticToEcl(0, Math.PI / 2);                   /* north galactic pole */
  const dd = ex.x * ez.x + ex.y * ez.y + ex.z * ez.z;
  ez = norm({ x: ez.x - dd * ex.x, y: ez.y - dd * ex.y, z: ez.z - dd * ex.z });
  const ey = cross(ez, ex);                                 /* right-handed with ez up */
  const toAU = function (u, v, w) {                         /* kpc in the galaxy frame -> AU, ecliptic */
    return { x: c.x + (u * ex.x + v * ey.x + w * ez.x) * KPC,
             y: c.y + (u * ex.y + v * ey.y + w * ez.y) * KPC,
             z: c.z + (u * ex.z + v * ey.z + w * ez.z) * KPC };
  };
  const rng = mulberry32(20260925);
  const g1 = function () { return Math.sqrt(-2 * Math.log(1 - rng() + 1e-12)) * Math.cos(TAU * rng()); };
  const BAR = -30 * DEG, cb = Math.cos(BAR), sb = Math.sin(BAR);

  /* arms: r = r0 · exp(k·θ), θ counter-clockwise from the Sun's azimuth (radians) */
  const ARMS = [
    { name: 'Scutum-Centaurus', r0: 4.5,  k: Math.tan(12 * DEG), t0: -30 * DEG,  t1: 370 * DEG, w: 0.40, amp: 1.0,  hot: 0.55, n: 1900 },
    { name: 'Perseus',          r0: 9.9,  k: Math.tan(10 * DEG), t0: -210 * DEG, t1: 150 * DEG, w: 0.42, amp: 1.0,  hot: 0.55, n: 1900 },
    { name: 'Sagittarius-Carina', r0: 6.4, k: Math.tan(12 * DEG), t0: -120 * DEG, t1: 220 * DEG, w: 0.30, amp: 0.7, hot: 0.45, n: 1200 },
    { name: 'Norma-Outer',      r0: 12.7, k: Math.tan(12 * DEG), t0: -300 * DEG, t1: 40 * DEG,  w: 0.30, amp: 0.55, hot: 0.40, n: 900 },
    { name: 'Local (Orion)',    r0: 8.35, k: Math.tan(11 * DEG), t0: -28 * DEG,  t1: 42 * DEG,  w: 0.20, amp: 0.45, hot: 0.50, n: 380, spur: true }
  ];
  const armFade = function (arm, th, r) {
    const f1 = arm.spur ? 1 : clamp((th - arm.t0) / (25 * DEG), 0, 1);   /* emerges from the bar */
    return f1 * clamp((15.5 - r) / 3.5, 0, 1);                           /* dies at the disc edge */
  };

  /* ---- particles ---- */
  const parts = [];
  const addP = function (u, v, w, cls, al) { parts.push({ u: u, v: v, w: w, c: cls, a: al }); };
  ARMS.forEach(function (arm) {
    for (let i = 0; i < arm.n * 1.6; i++) {
      const th = arm.t0 + rng() * (arm.t1 - arm.t0);
      const r = arm.r0 * Math.exp(arm.k * th);
      if (r > 15.5 || r < 3.2) continue;
      const f = armFade(arm, th, r);
      if (rng() > f) continue;
      const rr = r + g1() * arm.w, z = g1() * 0.09;
      const t = rng();
      addP(rr * Math.cos(th), rr * Math.sin(th), z, t < arm.hot ? 0 : (t < arm.hot + 0.08 ? 3 : 1), 0.55 + 0.45 * rng());
    }
  });
  for (let i = 0; i < 2800; i++) {                          /* the smooth disc, scale length 2.6 kpc */
    const r = -2.6 * (Math.log(1 - rng()) + Math.log(1 - rng()));
    if (r > 15 || r < 1.2) continue;
    const th = rng() * TAU;
    addP(r * Math.cos(th), r * Math.sin(th), g1() * (0.26 + 0.03 * r), rng() < 0.55 ? 2 : 1, 0.3 + 0.4 * rng());
  }
  for (let i = 0; i < 1300; i++) {                          /* bar and boxy bulge */
    let u, v, w;
    if (rng() < 0.55) { u = clamp(g1() * 1.6, -3.6, 3.6); v = g1() * 0.5; w = g1() * 0.35; }
    else { const rr = Math.abs(g1()) * 0.9, ph = rng() * TAU; u = rr * Math.cos(ph); v = rr * Math.sin(ph) * 0.85; w = g1() * 0.55; }
    addP(u * cb - v * sb, u * sb + v * cb, w, 2, 0.45 + 0.5 * rng());
  }
  for (let i = 0; i < 260; i++) {                           /* a sparse stellar halo */
    const r = 2.5 + Math.abs(g1()) * 7, ph = rng() * TAU, ct = 2 * rng() - 1, st = Math.sqrt(1 - ct * ct);
    addP(r * st * Math.cos(ph), r * st * Math.sin(ph), r * ct, 2, 0.22);
  }
  /* pack, sorted by colour bucket so fillStyle changes only a few times per frame */
  parts.forEach(function (p) { p.b = p.c * 2 + (p.a < 0.6 ? 0 : 1); });
  parts.sort(function (p, q) { return p.b - q.b; });
  const n = parts.length;
  const px = new Float32Array(n), py = new Float32Array(n), pz = new Float32Array(n);
  const bucketStart = new Int32Array(9);
  let cur = -1;
  parts.forEach(function (p, i) {
    const P = toAU(p.u, p.v, p.w); px[i] = P.x; py[i] = P.y; pz[i] = P.z;
    while (cur < p.b) { cur++; bucketStart[cur] = i; }
  });
  while (cur < 8) { cur++; bucketStart[cur] = n; }
  const COLS = ['170,190,255', '236,236,250', '255,224,176', '255,170,190'];
  const styles = [];
  for (let b = 0; b < 8; b++) styles.push({ col: COLS[b >> 1], a: (b & 1) ? 0.95 : 0.55 });

  /* ---- soft light and dust ---- */
  const blobs = [];
  const addB = function (u, v, w, rKpc, col, al, dust) {
    const P = toAU(u, v, w);
    blobs.push({ x: P.x, y: P.y, z: P.z, r: rKpc * KPC, col: col, a: al, dust: !!dust });
  };
  addB(0, 0, 0, 0.5, '255,240,215', 0.32);
  for (let i = 0; i < 22; i++) {
    const u = g1() * 1.2, v = g1() * 0.8, w = g1() * 0.5;
    addB(u * cb - v * sb, u * sb + v * cb, w, 0.9 + rng() * 0.6, '255,222,170', 0.07);
  }
  for (let i = 0; i < 10; i++) { const s = -3.2 + 6.4 * i / 9; addB(s * cb, s * sb, 0, 1.0, '255,215,160', 0.055); }
  ARMS.forEach(function (arm) {
    const wf = arm.w / 0.4;
    for (let th = arm.t0; th < arm.t1; th += 22 * DEG) {
      const r = arm.r0 * Math.exp(arm.k * th);
      const f = armFade(arm, th, r); if (r > 14.5 || f < 0.15) continue;
      const rr = r + g1() * 0.15;
      addB(rr * Math.cos(th), rr * Math.sin(th), g1() * 0.05, (0.9 + rng() * 0.5) * wf, arm.spur ? '210,215,255' : '190,205,255', 0.075 * arm.amp * f);
    }
    for (let th = arm.t0; th < arm.t1; th += 13 * DEG) {    /* dust hugs the inner edge of each arm */
      const r = arm.r0 * Math.exp(arm.k * th);
      const f = armFade(arm, th, r); if (r > 13 || f < 0.2) continue;
      const rd = r * 0.92 - 0.1 + g1() * 0.08;
      addB(rd * Math.cos(th), rd * Math.sin(th), g1() * 0.04, (0.55 + rng() * 0.3) * wf, '14,9,6', 0.40 * arm.amp * f, true);
    }
  });
  for (let i = 0; i < 26; i++) {
    const r = -3.2 * (Math.log(1 - rng()) + Math.log(1 - rng()));
    if (r > 13 || r < 2) { i--; continue; }
    const th = rng() * TAU;
    addB(r * Math.cos(th), r * Math.sin(th), g1() * 0.15, 2.2 + rng() * 1.2, '232,222,205', 0.04);
  }
  for (let i = 0; i < 22; i++) {                            /* dust in the midplane */
    const r = 2.5 + rng() * 7, th = rng() * TAU;
    addB(r * Math.cos(th), r * Math.sin(th), 0, 1.0 + rng() * 0.6, '14,9,6', 0.20, true);
  }
  [[8, 0.018], [11, 0.014], [14, 0.010], [17, 0.007], [20, 0.005]].forEach(function (h) { addB(0, 0, 0, h[0], '205,205,235', h[1]); });
  const order = blobs.map(function (_, i) { return i; });
  const dist2 = new Float64Array(blobs.length);
  /* where to write each arm's name when the galaxy is seen from outside */
  const armLabels = [];
  const armAt = function (i, thDeg) { const arm = ARMS[i], th = thDeg * DEG, r = arm.r0 * Math.exp(arm.k * th); return toAU(r * Math.cos(th), r * Math.sin(th), 0); };
  armLabels.push({ arm: 0, pos: armAt(0, 120) }, { arm: 1, pos: armAt(1, 40) }, { arm: 2, pos: armAt(2, -60) },
                 { arm: 3, pos: armAt(3, -240) }, { arm: 3, pos: armAt(3, 20), alt: 'Outer Arm' }, { arm: 4, pos: armAt(4, 32) },
                 { arm: null, pos: toAU(3.0 * cb, 3.0 * sb, 0), alt: 'Galactic bar' });

  /* the selectable object: lives at the centre, framed from above the plane on the Sun's side */
  const view = norm({ x: ex.x * Math.cos(38 * DEG) + ez.x * Math.sin(38 * DEG),
                      y: ex.y * Math.cos(38 * DEG) + ez.y * Math.sin(38 * DEG),
                      z: ex.z * Math.cos(38 * DEG) + ez.z * Math.sin(38 * DEG) });
  const radius = 50000 * LY_AU;
  const o = addObject({
    id: MILKY_WAY.id, name: MILKY_WAY.name, type: MILKY_WAY.type, kind: 'region', category: 'galaxies',
    color: MILKY_WAY.color, radiusAU: radius, regionRadius: radius,
    fixed: { x: c.x, y: c.y, z: c.z }, data: MILKY_WAY, prio: 90, minPx: 0, noRender: true,
    centre: true, viewFill: 0.75,
    view: { yaw: Math.atan2(view.y, view.x), pitch: Math.asin(clamp(view.z, -1, 1)) }
  });
  o.pos = o.fixed;
  return { centre: c, ex: ex, ey: ey, ez: ez, radius: radius, KPC: KPC, n: n, px: px, py: py, pz: pz, arms: ARMS, armLabels: armLabels, toAU: toAU,
           bucketStart: bucketStart, styles: styles, blobs: blobs, order: order, dist2: dist2, obj: o };
})();


/* --- The cosmic address: arms, the Local Group, the Virgo Supercluster,
   Laniakea and the observable universe as selectable regions --- */
const STRUCT = (function () {
  if (!MW || typeof STRUCTURE === 'undefined') return null;
  const MLY = 1e6 * LY_AU;
  const regions = [];
  const view38 = function (n) {                          /* a viewpoint 38° above a plane with normal n, on our side */
    const toUs = norm({ x: -n.cx, y: -n.cy, z: -n.cz });
    const nn = n.normal;
    const dd = toUs.x * nn.x + toUs.y * nn.y + toUs.z * nn.z;
    const u = norm({ x: toUs.x - dd * nn.x, y: toUs.y - dd * nn.y, z: toUs.z - dd * nn.z });
    const v = norm({ x: u.x * Math.cos(38 * DEG) + nn.x * Math.sin(38 * DEG), y: u.y * Math.cos(38 * DEG) + nn.y * Math.sin(38 * DEG), z: u.z * Math.cos(38 * DEG) + nn.z * Math.sin(38 * DEG) });
    return { yaw: Math.atan2(v.y, v.x), pitch: Math.asin(clamp(v.z, -1, 1)) };
  };
  const mkRegion = function (cfg, centre, radiusAU, normal, level) {
    const o = addObject({
      id: cfg.id, name: cfg.name, type: 'Large-scale structure', kind: 'region', category: 'structure',
      color: cfg.color, radiusAU: radiusAU, regionRadius: radiusAU,
      fixed: { x: centre.x, y: centre.y, z: centre.z }, data: cfg, prio: 85, minPx: 0, noRender: true,
      centre: true, viewFill: 0.75, structLevel: level,
      view: view38({ cx: centre.x, cy: centre.y, cz: centre.z, normal: normal })
    });
    o.pos = o.fixed;
    regions.push({ obj: o, centre: o.pos, R: radiusAU, normal: normal, col: cfg.color });
    return o;
  };
  /* arms are places within the galaxy; selecting one shows the galaxy with that arm named */
  STRUCTURE.arms.forEach(function (cfg) {
    const o = addObject({
      id: cfg.id, name: cfg.name, type: cfg.arm == null ? 'Galactic centre' : 'Spiral arm', kind: 'region', category: 'structure',
      color: cfg.color, radiusAU: MW.radius, regionRadius: MW.radius,
      fixed: { x: MW.centre.x, y: MW.centre.y, z: MW.centre.z }, data: cfg, prio: 80, minPx: 0, noRender: true,
      centre: true, viewFill: 0.75, structLevel: 0, armIndex: cfg.arm,
      view: MW.obj.view
    });
    o.pos = o.fixed;
  });
  MW.obj.structLevel = 1;
  const m31 = byId.m31, m87 = byId.m87;
  /* Local Group: centred between the two big spirals */
  const lgC = m31 ? { x: (MW.centre.x + m31.pos.x) / 2, y: (MW.centre.y + m31.pos.y) / 2, z: (MW.centre.z + m31.pos.z) / 2 } : MW.centre;
  mkRegion(STRUCTURE.localgroup, lgC, STRUCTURE.localgroup.radiusMly * MLY, MW.ez, 2);
  /* Virgo (Local) Supercluster: a disc in the supergalactic plane, centred on the Virgo Cluster */
  const sgNorth = galacticToEcl(47.37 * DEG, 6.32 * DEG);
  const vC = m87 ? m87.pos : norm(raDecToVec(12.5137, 12.391, 1));
  const vCentre = m87 ? { x: vC.x, y: vC.y, z: vC.z } : { x: vC.x * 54 * MLY, y: vC.y * 54 * MLY, z: vC.z * 54 * MLY };
  const vsc = mkRegion(STRUCTURE.virgosc, vCentre, STRUCTURE.virgosc.radiusMly * MLY, sgNorth, 3);
  regions[regions.length - 1].disc = true;
  /* Laniakea: a basin of attraction around the Great Attractor (Norma cluster direction) */
  const ga = galacticToEcl(325.3 * DEG, -7.3 * DEG);
  mkRegion(STRUCTURE.laniakea, { x: ga.x * 250 * MLY, y: ga.y * 250 * MLY, z: ga.z * 250 * MLY }, STRUCTURE.laniakea.radiusMly * MLY, sgNorth, 4);
  /* the observable universe: centred on the observer */
  mkRegion(STRUCTURE.observable, { x: 0, y: 0, z: 0 }, STRUCTURE.observable.radiusMly * MLY, MW.ez, 5);
  return { regions: regions, chain: ['milkyway', 'localgroup', 'virgosc', 'laniakea', 'observable'], lgC: lgC, vCentre: vCentre };
})();

/* Which arm (or what else) a point in the galaxy belongs to */
function galaxyPlace(pos) {
  const K = MW.KPC, du = pos.x - MW.centre.x, dv = pos.y - MW.centre.y, dw = pos.z - MW.centre.z;
  const u = (du * MW.ex.x + dv * MW.ex.y + dw * MW.ex.z) / K;
  const v = (du * MW.ey.x + dv * MW.ey.y + dw * MW.ey.z) / K;
  const w = (du * MW.ez.x + dv * MW.ez.y + dw * MW.ez.z) / K;
  const r = Math.hypot(u, v), th = Math.atan2(v, u);
  if (r < 3.3) return { id: 'bar', name: 'Galactic bar and bulge' };
  if (Math.abs(w) > 2.5) return { id: null, name: 'Galactic halo' };
  const ids = ['arm-scutum', 'arm-perseus', 'arm-sagittarius', 'arm-norma', 'arm-orion'];
  const hits = [];
  MW.arms.forEach(function (arm, i) {
    for (let m = -2; m <= 2; m++) {
      const t = th + m * TAU;
      if (t < arm.t0 || t > arm.t1) continue;
      const rr = arm.r0 * Math.exp(arm.k * t);
      hits.push({ i: i, gap: Math.abs(r - rr) * Math.cos(Math.atan(arm.k)), inner: rr < r });
    }
  });
  hits.sort(function (p, q) { return p.gap - q.gap; });
  if (!hits.length) return { id: null, name: 'Outer disc' };
  if (hits[0].gap < 0.7) return { id: ids[hits[0].i], name: byId[ids[hits[0].i]].name };
  const inner = hits.find(function (h) { return h.inner; }), outer = hits.find(function (h) { return !h.inner; });
  if (inner && outer) return { id: null, name: 'Between the ' + byId[ids[inner.i]].name + ' and the ' + byId[ids[outer.i]].name };
  return { id: null, name: 'Near the ' + byId[ids[hits[0].i]].name };
}

/* The chain of structures an object belongs to, innermost first */
function cosmicAddress(o) {
  if (!MW || !STRUCT) return [];
  const steps = [];
  const push = function (id, name) { if (id && !byId[id]) return; steps.push({ id: id, name: name || byId[id].name }); };
  const k = o.kind;
  if (k === 'constellation') return steps;                 /* a pattern on the sky, not a place */
  if (k === 'deepfield') { push('observable', 'Observable universe, along this line of sight'); return steps; }
  const MLY = 1e6 * LY_AU;
  const tail = function (pos) {
    if (o.id !== 'localgroup' && o.structLevel == null || (o.structLevel != null && o.structLevel < 2)) {
      const dLG = Math.hypot(pos.x - STRUCT.lgC.x, pos.y - STRUCT.lgC.y, pos.z - STRUCT.lgC.z);
      if (dLG < 5 * MLY) push('localgroup');
    }
    if (o.structLevel == null || o.structLevel < 3) {
      const dV = Math.hypot(pos.x - STRUCT.vCentre.x, pos.y - STRUCT.vCentre.y, pos.z - STRUCT.vCentre.z);
      if (dV < 60 * MLY) push('virgosc');
    }
    if (o.structLevel == null || o.structLevel < 4) {
      const L = byId.laniakea, dL = Math.hypot(pos.x - L.pos.x, pos.y - L.pos.y, pos.z - L.pos.z);
      if (dL < 262 * MLY) push('laniakea'); else steps.push({ id: null, name: 'Beyond Laniakea' });
    }
    if (o.structLevel == null || o.structLevel < 5) push('observable');
  };
  if (o.structLevel != null) {                              /* a structure: everything above it */
    if (o.structLevel === 0) push('milkyway');
    tail(o.structLevel <= 1 ? MW.centre : o.pos);
    return steps;
  }
  const solar = k === 'sun' || k === 'planet' || k === 'moon' || k === 'dwarf' || k === 'probe' || o.category === 'satellites' || o.category === 'belts';
  if (solar) {
    if (o.category === 'satellites' || (k === 'probe' && (o.mode === 'earthorbit' || o.mode === 'gp'))) push('earth');
    else if (o.parent && byId[o.parent] && k !== 'exoplanet') push(o.parent);
    if (k !== 'sun') push('sun', 'Solar System');
    push('arm-orion'); push('milkyway'); tail(MW.centre); return steps;
  }
  if (k === 'ngc' || (o.dir && !o.pos)) { push('milkyway', 'Milky Way, distance not known'); tail(MW.centre); return steps; }
  let pos = o.pos;
  if (k === 'exoplanet' && o.parent && byId[o.parent]) { push(o.parent); pos = byId[o.parent].pos; }
  if (!pos) return steps;
  const isGalaxy = o.category === 'galaxies' || (o.sub === 'galaxy' || o.sub === 'spiral');
  const dGC = Math.hypot(pos.x - MW.centre.x, pos.y - MW.centre.y, pos.z - MW.centre.z);
  if (!isGalaxy && dGC < 70000 * LY_AU) {
    const pl = galaxyPlace(pos);
    steps.push(pl);
    push('milkyway'); tail(MW.centre); return steps;
  }
  if (o.data && o.data.within && o.data.within !== 'laniakea' && o.data.within !== 'virgosc') push(o.data.within);
  tail(pos);
  return steps;
}

/* Outlines of the big structures, once the camera is outside them */
function drawWall(g, dist) {
  const sel = selectedId === g.obj.id;
  /* a wall only makes sense from outside it: fade in once the camera is well beyond the local universe */
  const camR = Math.hypot(camPos.x, camPos.y, camPos.z), dW = Math.hypot(g.centre.x, g.centre.y, g.centre.z);
  let al = sel ? 1 : clamp((camR / dW - 0.45) / 0.5, 0, 1);
  if (dist < g.R * 0.4) al = 0;
  if (al <= 0.02) return;
  const pts = [];
  for (let k = 0; k < g.path.length; k++) { const q = project(g.path[k]); if (!q) return; pts.push(q); }
  let len = 0; for (let k = 1; k < pts.length; k++) len += Math.hypot(pts[k].x - pts[k - 1].x, pts[k].y - pts[k - 1].y);
  if (len < 50) return;
  al *= clamp((len - 50) / 120, 0, 1) * clamp(4 * Math.max(W, H) / len, 0, 1);
  if (al <= 0.02) return;
  const rgb = hexA(g.col, 1).replace('rgba(', '').replace(',1)', '');
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  const C = g.cloud, fx = fwd.x, fy = fwd.y, fz = fwd.z;
  for (let k = 0; k < C.n; k++) {
    const dx = C.x[k] - camPos.x, dy = C.y[k] - camPos.y, dz = C.z[k] - camPos.z;
    const vz = dx * fx + dy * fy + dz * fz; if (vz <= 0) continue;
    const sx = cx + (dx * right.x + dy * right.y + dz * right.z) * focal / vz;
    const sy = cy - (dx * up.x + dy * up.y + dz * up.z) * focal / vz;
    if (sx < -5 || sx > W + 5 || sy < -5 || sy > H + 5) continue;
    const b = C.b[k];
    ctx.fillStyle = 'rgba(' + rgb + ',' + ((0.25 + 0.6 * b) * al).toFixed(3) + ')';
    const s = 0.9 + b * 1.6;
    ctx.fillRect(sx - s / 2, sy - s / 2, s, s);
  }
  ctx.globalCompositeOperation = 'source-over';
  if (sel) {
    ctx.setLineDash([5, 6]); ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(' + rgb + ',0.6)';
    ctx.beginPath(); ctx.moveTo(pts[0].x, pts[0].y);
    for (let k = 1; k < pts.length; k++) ctx.lineTo(pts[k].x, pts[k].y);
    ctx.stroke(); ctx.setLineDash([]);
  }
  const m = pts[Math.floor(pts.length / 2)];
  const txt = g.obj.name;
  const tw = ctx.measureText(txt).width;
  const box = { x: m.x - tw / 2 - 4, y: m.y - 30, w: tw + 8, h: 16 };
  if (!overlapsLabel(box) || sel) {
    backdropLabelBoxes.push(box);
    ctx.fillStyle = 'rgba(' + rgb + ',' + (al * (sel ? 1 : 0.85)).toFixed(3) + ')';
    ctx.fillText(txt, m.x, m.y - 22);
  }
  ctx.restore();
}
function drawRegions() {
  if (!STRUCT) return;
  const maxWH = Math.max(W, H);
  const camR = Math.hypot(camPos.x, camPos.y, camPos.z);
  ctx.save();
  ctx.font = '600 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  if ('letterSpacing' in ctx) ctx.letterSpacing = '1px';
  ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  for (let i = 0; i < STRUCT.regions.length; i++) {
    const g = STRUCT.regions[i], o = g.obj;
    if (deep.on && ((deep.shot !== 'cosmos' && deep.shot !== 'lg') || (o.id !== 'localgroup' && o.id !== 'observable' && deep.horizonFade < 0.5))) continue;
    const dx = g.centre.x - camPos.x, dy = g.centre.y - camPos.y, dz = g.centre.z - camPos.z;
    const dist = Math.hypot(dx, dy, dz);
    if (g.path) { drawWall(g, dist); continue; }
    if (dist <= g.R * 1.03) continue;                       /* inside it */
    const c = project(g.centre); if (!c) continue;
    const sr = (g.R / Math.sqrt(dist * dist - g.R * g.R)) * focal;
    if (sr < 7 || sr > 4 * maxWH) continue;
    const sel = selectedId === o.id;
    let al = clamp((sr - 7) / 25, 0, 1);
    if (sr > 1.5 * maxWH) al *= clamp((4 * maxWH - sr) / (2.5 * maxWH), 0, 1);
    if (al <= 0.01) continue;
    const rgb = hexA(g.col, 1).replace('rgba(', '').replace(',1)', '');
    if (opts.structures === 'off' && !sel) continue;
    const shells = opts.structures !== 'outlines';
    ctx.setLineDash(sel || shells ? [] : [6, 5]);
    ctx.lineWidth = sel ? 1.6 : 1;
    ctx.strokeStyle = 'rgba(' + rgb + ',' + (al * (sel ? 0.95 : (shells ? 0.18 : 0.55))).toFixed(3) + ')';
    if (g.disc) {
      /* a flattened structure: draw its rim in 3D */
      const n = g.normal, t = Math.abs(n.z) < 0.9 ? { x: 0, y: 0, z: 1 } : { x: 1, y: 0, z: 0 };
      const u1 = norm(cross(n, t)), u2 = cross(n, u1);
      ctx.beginPath(); let ok = true;
      for (let k = 0; k <= 72; k++) {
        const ang = k / 72 * TAU;
        const q = project({ x: g.centre.x + g.R * (Math.cos(ang) * u1.x + Math.sin(ang) * u2.x),
                            y: g.centre.y + g.R * (Math.cos(ang) * u1.y + Math.sin(ang) * u2.y),
                            z: g.centre.z + g.R * (Math.cos(ang) * u1.z + Math.sin(ang) * u2.z) });
        if (!q) { ok = false; break; }
        if (k === 0) ctx.moveTo(q.x, q.y); else ctx.lineTo(q.x, q.y);
      }
      if (ok) {
        if (shells) {                                      /* a flattened shell: a soft rim and a faint film */
          ctx.save(); ctx.globalCompositeOperation = 'lighter'; ctx.lineJoin = 'round';
          [[18, 0.05], [8, 0.10], [2.5, 0.28]].forEach(function (s) { ctx.lineWidth = s[0]; ctx.strokeStyle = 'rgba(' + rgb + ',' + (s[1] * al).toFixed(3) + ')'; ctx.stroke(); });
          ctx.fillStyle = 'rgba(' + rgb + ',' + (0.045 * al).toFixed(3) + ')'; ctx.fill();
          ctx.restore();
        }
        if (!shells || sel) ctx.stroke();
      }
    } else {
      if (shells) {                                        /* a bubble: bright at the rim, clear in the middle */
        const gr = ctx.createRadialGradient(c.x, c.y, sr * 0.72, c.x, c.y, sr * 1.03);
        const a0 = g.isVoid ? 0.5 : 1;
        gr.addColorStop(0, 'rgba(' + rgb + ',0)'); gr.addColorStop(0.72, 'rgba(' + rgb + ',' + (0.07 * al * a0).toFixed(3) + ')');
        gr.addColorStop(0.95, 'rgba(' + rgb + ',' + (0.34 * al * a0).toFixed(3) + ')'); gr.addColorStop(0.985, 'rgba(' + rgb + ',' + (0.42 * al * a0).toFixed(3) + ')'); gr.addColorStop(1, 'rgba(' + rgb + ',0)');
        ctx.save(); ctx.globalCompositeOperation = 'lighter';
        ctx.fillStyle = gr; ctx.beginPath(); ctx.arc(c.x, c.y, sr * 1.03, 0, TAU); ctx.fill();
        ctx.restore();
        if (sel) { ctx.beginPath(); ctx.arc(c.x, c.y, sr, 0, TAU); ctx.stroke(); }
      } else { ctx.beginPath(); ctx.arc(c.x, c.y, sr, 0, TAU); ctx.stroke(); }
      if (g.isVoid) {                                       /* a void: darker than its surroundings */
        const gr = ctx.createRadialGradient(c.x, c.y, 0, c.x, c.y, sr);
        gr.addColorStop(0, 'rgba(0,0,0,' + (0.45 * al).toFixed(3) + ')'); gr.addColorStop(1, 'rgba(0,0,0,0)');
        ctx.fillStyle = gr; ctx.beginPath(); ctx.arc(c.x, c.y, sr, 0, TAU); ctx.fill();
      }
      if (o.id === 'observable') {                          /* the CMB: a faint glow just inside the edge */
        const gr = ctx.createRadialGradient(c.x, c.y, sr * 0.9, c.x, c.y, sr);
        gr.addColorStop(0, 'rgba(' + rgb + ',0)'); gr.addColorStop(1, 'rgba(' + rgb + ',' + (0.16 * al).toFixed(3) + ')');
        ctx.fillStyle = gr; ctx.beginPath(); ctx.arc(c.x, c.y, sr, 0, TAU); ctx.fill();
      }
    }
    ctx.setLineDash([]);
    const ly = clamp(c.y - sr - 9, 26, H - 26);
    const txt = o.name + ' · ' + fmtDist(g.R * 2) + ' across';
    const tw = ctx.measureText(txt).width;
    const box = { x: c.x - tw / 2 - 4, y: ly - 8, w: tw + 8, h: 16 };
    if (c.x - tw / 2 > 10 && c.x + tw / 2 < W - 10 && (!overlapsLabel(box) || sel)) {
      backdropLabelBoxes.push(box);
      ctx.fillStyle = 'rgba(' + rgb + ',' + (al * (sel ? 1 : 0.8)).toFixed(3) + ')';
      ctx.fillText(txt, c.x, ly);
    }
  }
  /* far out, the catalogue is a speck: say what it is */
  if (camR > 4e9 * LY_AU) {
    const c0 = project({ x: 0, y: 0, z: 0 });
    if (c0 && bulkGal) {
      ctx.textAlign = 'left';
      ctx.font = '500 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
      ctx.fillStyle = 'rgba(220,215,245,0.6)';
      ctx.fillText('the ' + bulkGal.n.toLocaleString() + ' mapped galaxies', c0.x + 12, c0.y);
    }
  }
  ctx.restore();
}

/* --- Background sky (directions only, effectively at infinity) --- */
const skyStars = (function () {
  const rng = mulberry32(20260902);
  const n = 3600, arr = [];
  for (let k = 0; k < n; k++) {
    let dir, mag;
    if (k < n * 0.55) {                       // Milky Way band
      const l = rng() * TAU;
      const b = (rng() + rng() + rng() - 1.5) * 11 * DEG;
      dir = galacticToEcl(l, b);
      mag = 0.28 + rng() * 0.72;
    } else {                                  // isotropic field
      const u = rng() * 2 - 1, ph = rng() * TAU, s = Math.sqrt(1 - u * u);
      dir = { x: s * Math.cos(ph), y: s * Math.sin(ph), z: u };
      mag = 0.18 + rng() * 0.6;
    }
    arr.push({ d: dir, m: mag, s: rng() < 0.06 ? 2 : 1, tint: rng() });
  }
  arr.sort((a, b) => a.m - b.m);
  return arr;
})();

/* The Milky Way: a volumetric band built from many soft blobs placed in
   galactic coordinates, with dark dust-lane blobs painted over the glow.
   Sizes are angular, so the band scales correctly with the window. */
const mwBlobs = (function () {
  const rng = mulberry32(77713);
  const arr = [];
  function push(l, b, angDeg, a, col, dust) {
    arr.push({ d: galacticToEcl(l, b), ang: angDeg * DEG, a: a, col: col, dust: dust });
  }
  /* glowing band — thicker and warmer toward the galactic centre */
  for (let k = 0; k < 430; k++) {
    const l = rng() * TAU;
    const dl = Math.min(l, TAU - l);
    const centre = Math.pow(Math.max(0, 1 - dl / Math.PI), 1.8);
    const sigma = (2.5 + centre * 7.5) * DEG;
    const b = (rng() + rng() + rng() - 1.5) * sigma;
    const ang = 1.6 + rng() * 5.5 + centre * 7;
    const a = (0.045 + rng() * 0.060) * (0.30 + centre * 2.0);
    const col = centre > 0.55
      ? 'rgb(' + (238 - (rng() * 18 | 0)) + ',' + (216 - (rng() * 26 | 0)) + ',' + (184 - (rng() * 32 | 0)) + ')'
      : 'rgb(' + (194 - (rng() * 22 | 0)) + ',' + (208 - (rng() * 18 | 0)) + ',' + (236 - (rng() * 14 | 0)) + ')';
    push(l, b, ang, a, col, false);
  }
  /* the bulge itself: a few large, soft, warm masses */
  for (let k = 0; k < 30; k++) {
    const l = (rng() - 0.5) * 48 * DEG;
    const b = (rng() + rng() - 1) * 9 * DEG;
    push(l, b, 9 + rng() * 16, 0.085 + rng() * 0.075, 'rgb(243,221,188)', false);
  }
  /* dark rift: dust lanes hugging the plane, painted over the glow */
  for (let k = 0; k < 230; k++) {
    const l = rng() * TAU;
    const dl = Math.min(l, TAU - l);
    const centre = Math.pow(Math.max(0, 1 - dl / Math.PI), 1.5);
    const b = (rng() + rng() - 1) * 1.6 * DEG - 0.3 * DEG;
    const ang = 1.2 + rng() * 4.4 + centre * 3.5;
    push(l, b, ang, (0.10 + rng() * 0.20) * (0.4 + centre), 'rgb(11,8,7)', true);
  }
  return arr;
})();

/* colour helper: '#rrggbb' + alpha -> rgba() */
function hexA(hex, a) {
  return 'rgba(' + parseInt(hex.slice(1, 3), 16) + ',' + parseInt(hex.slice(3, 5), 16) +
         ',' + parseInt(hex.slice(5, 7), 16) + ',' + a + ')';
}

/* ============================================================
   3b. Bundled catalogues
   Packed binary in .js files rather than JSON, because file:// blocks
   fetch() and XHR — a classic <script> tag is the only thing that loads
   from disk without a server.
   ============================================================ */
const PC_AU = 206264.806;

function b64buf(str) {
  const bin = atob(str), n = bin.length, u8 = new Uint8Array(n);
  for (let i = 0; i < n; i++) u8[i] = bin.charCodeAt(i);
  return u8.buffer;
}

let bulkStars = null, bulkDso = null, bulkExo = null;

if (typeof HYG_BIN !== 'undefined') {
  const b = b64buf(HYG_BIN), n = HYG_N;
  const sx = new Float32Array(b, 0, n), sy = new Float32Array(b, 4 * n, n), sz = new Float32Array(b, 8 * n, n);
  for (let i = 0; i < n; i++) { sx[i] *= PC_AU; sy[i] *= PC_AU; sz[i] *= PC_AU; }   /* parsecs -> AU */
  bulkStars = { n: n, x: sx, y: sy, z: sz,
                mag: new Int16Array(b, 12 * n, n), ci: new Int16Array(b, 14 * n, n), named: {} };
  HYG_NAMES.forEach(function (e) { bulkStars.named[e[0]] = e; });
  /* default cut: magnitude 7.5 — a dark-sky naked-eye view with a little extra depth.
     The array is sorted by magnitude, so the cut is just an index. */
  { let lo = 0, hi = n; while (lo < hi) { const m = (lo + hi) >> 1; if (bulkStars.mag[m] <= 750) lo = m + 1; else hi = m; } bulkStars.defaultN = lo; }
}

if (typeof NGC_BIN !== 'undefined') {
  const b = b64buf(NGC_BIN), n = NGC_N;
  bulkDso = { n: n,
    x: new Int16Array(b, 0, n), y: new Int16Array(b, 2 * n, n), z: new Int16Array(b, 4 * n, n),
    maj: new Uint16Array(b, 6 * n, n), mag: new Int16Array(b, 8 * n, n),
    type: new Uint8Array(b, 10 * n, n),
    names: NGC_NAMES.split('|'), common: NGC_COMMON.split('|'), typeNames: NGC_TYPES };
  bulkDso.lname = bulkDso.names.map(function (v) { return v.toLowerCase(); });
  bulkDso.lcommon = bulkDso.common.map(function (v) { return v.toLowerCase(); });
}

if (typeof EXO_BIN !== 'undefined') {
  const b = b64buf(EXO_BIN), n = EXO_N;
  bulkExo = { n: n,
    dist: new Uint32Array(b, 0, n), a: new Uint32Array(b, 4 * n, n), per: new Uint32Array(b, 8 * n, n),
    x: new Int16Array(b, 12 * n, n), y: new Int16Array(b, 14 * n, n), z: new Int16Array(b, 16 * n, n),
    rad: new Uint16Array(b, 18 * n, n), year: new Uint16Array(b, 20 * n, n),
    names: EXO_NAMES.split('|'), hosts: EXO_HOSTS.split('|') };
  bulkExo.lname = bulkExo.names.map(function (v) { return v.toLowerCase(); });
}


/* ---- JPL Small-Body Database: asteroids and comets ----
   Layout (32-bit first, then 16-bit, then 8-bit, so every view is aligned):
     f32 a[n] | u16 e*65535 | u16 i | u16 om | u16 w | u16 M0 (angles as 0..65535 = 0..360)
     u16 diam*10 km (0 = unknown) | i16 H*100 | u8 class
   All elements are propagated to one common epoch SBDB_EPOCH at build time. */
let bulkSbdb = null;
if (typeof SBDB_BIN !== 'undefined') {
  const b = b64buf(SBDB_BIN), n = SBDB_N;
  const A = new Float32Array(b, 0, n);
  const E = new Uint16Array(b, 4 * n, n), I = new Uint16Array(b, 6 * n, n), OM = new Uint16Array(b, 8 * n, n);
  const Wp = new Uint16Array(b, 10 * n, n), M0 = new Uint16Array(b, 12 * n, n), DI = new Uint16Array(b, 14 * n, n);
  const Hm = new Int16Array(b, 16 * n, n), CL = new Uint8Array(b, 18 * n, n);
  const K = 0.01720209895;                        /* Gaussian gravitational constant, rad/day */
  const ANG = TAU / 65536;
  /* precompute the orbital-plane basis vectors P (toward perihelion) and Q */
  const nrad = new Float32Array(n), ecc = new Float32Array(n), m0 = new Float32Array(n);
  const Px = new Float32Array(n), Py = new Float32Array(n), Pz = new Float32Array(n);
  const Qx = new Float32Array(n), Qy = new Float32Array(n), Qz = new Float32Array(n);
  for (let k = 0; k < n; k++) {
    const e = E[k] / 65535, inc = I[k] * ANG, O = OM[k] * ANG, w = Wp[k] * ANG;
    const cw = Math.cos(w), sw = Math.sin(w), cO = Math.cos(O), sO = Math.sin(O), cI = Math.cos(inc), sI = Math.sin(inc);
    Px[k] = cw * cO - sw * sO * cI;  Py[k] = cw * sO + sw * cO * cI;  Pz[k] = sw * sI;
    Qx[k] = -sw * cO - cw * sO * cI; Qy[k] = -sw * sO + cw * cO * cI; Qz[k] = cw * sI;
    ecc[k] = e; m0[k] = M0[k] * ANG;
    nrad[k] = K / Math.pow(A[k], 1.5);
  }
  bulkSbdb = { n: n, a: A, e: ecc, m0: m0, nrad: nrad, Px, Py, Pz, Qx, Qy, Qz, H: Hm, diam: DI, cls: CL,
               I: I, OM: OM, W: Wp,
               px: new Float32Array(n), py: new Float32Array(n), pz: new Float32Array(n),
               names: SBDB_NAMES.split('|'), keys: SBDB_KEYS.split('|'), epoch: SBDB_EPOCH, cursor: 0, lastJd: null };
  bulkSbdb.lname = bulkSbdb.names.map(function (v) { return v.toLowerCase(); });
  bulkSbdb.lkey = bulkSbdb.keys.map(function (v) { return v.toLowerCase(); });
  /* "important" = the large ones: H < 11 main belt, < 17 near-Earth, < 10.5 Trojans,
     < 8 trans-Neptunian, every comet. The rest are opt-in. */
  bulkSbdb.imp = new Uint8Array(n); bulkSbdb.nImp = 0;
  for (let k = 0; k < n; k++) {
    const c = bulkSbdb.cls[k], h = bulkSbdb.H[k];
    const ok = (c === 0 && h < 1100) || (c === 1 && h < 1700) || (c === 2 && h < 1050) || (c === 3 && h < 800) || c >= 4;
    if (ok) { bulkSbdb.imp[k] = 1; bulkSbdb.nImp++; }
  }
}

/* Two-body position for row k at Julian day jd, into the cache */
function sbdbSolve(S, k, jd) {
  const e = S.e[k];
  let M = S.m0[k] + S.nrad[k] * (jd - S.epoch);
  M -= TAU * Math.floor(M / TAU);
  let Ea = e < 0.8 ? M : Math.PI;
  for (let it = 0; it < 10; it++) {
    const d = (Ea - e * Math.sin(Ea) - M) / (1 - e * Math.cos(Ea));
    Ea -= d;
    if (Math.abs(d) < 1e-7) break;
  }
  const a = S.a[k];
  const xp = a * (Math.cos(Ea) - e), yp = a * Math.sqrt(1 - e * e) * Math.sin(Ea);
  S.px[k] = xp * S.Px[k] + yp * S.Qx[k];
  S.py[k] = xp * S.Py[k] + yp * S.Qy[k];
  S.pz[k] = xp * S.Pz[k] + yp * S.Qz[k];
}

/* Amortised: one sixth of the catalogue per frame, or everything after a time jump */
function sbdbUpdate(jd) {
  const S = bulkSbdb; if (!S) return;
  if (S.lastJd === null || Math.abs(jd - S.lastJd) > 20) {
    for (let k = 0; k < S.n; k++) sbdbSolve(S, k, jd);
    S.lastJd = jd; S.cursor = 0; return;
  }
  const step = Math.ceil(S.n / 6);
  const end = Math.min(S.n, S.cursor + step);
  for (let k = S.cursor; k < end; k++) sbdbSolve(S, k, jd);
  S.cursor = end >= S.n ? 0 : end;
  if (S.cursor === 0) S.lastJd = jd;
}

const SBDB_CLASS = [
  { key: 'mainbelt', label: 'Main-belt asteroid',   color: '#b9a98f', rgb: '185,169,143' },
  { key: 'neo',      label: 'Near-Earth asteroid',  color: '#ff8a5c', rgb: '255,138,92'  },
  { key: 'trojans',  label: 'Jupiter Trojan',       color: '#7fd0c0', rgb: '127,208,192' },
  { key: 'tno',      label: 'Trans-Neptunian object', color: '#8fa3e0', rgb: '143,163,224' },
  { key: 'comets',   label: 'Short-period comet',   color: '#9fe8ff', rgb: '159,232,255' },
  { key: 'comets',   label: 'Long-period comet',    color: '#c8f0ff', rgb: '200,240,255' }
];

function drawSmallBodies(sunScreen) {
  const S = bulkSbdb; if (!S || cam.dist > 6000 || deep.on) return;
  const px = camPos.x, py = camPos.y, pz = camPos.z;
  const fx = fwd.x, fy = fwd.y, fz = fwd.z, rx = right.x, ry = right.y, rz = right.z, ux = up.x, uy = up.y, uz = up.z;
  const on = [opts.mainbelt, opts.neo, opts.trojans, opts.tno, opts.comets, opts.comets];
  let cur = -1, drawn = 0;
  const budget = 45000;
  /* from far away a belt is a faint haze, not a clump of dots */
  const zoomFade = [clamp(9 / cam.dist, 0.18, 1), clamp(9 / cam.dist, 0.18, 1), clamp(14 / cam.dist, 0.18, 1),
                    clamp(110 / cam.dist, 0.18, 1), clamp(30 / cam.dist, 0.3, 1), clamp(30 / cam.dist, 0.3, 1)];
  ctx.save();
  for (let k = 0; k < S.n && drawn < budget; k++) {
    const c = S.cls[k]; if (!on[c]) continue;
    if (!opts.allSmallBodies && !S.imp[k]) continue;
    const dx = S.px[k] - px, dy = S.py[k] - py, dz = S.pz[k] - pz;
    const vz = dx * fx + dy * fy + dz * fz;
    if (vz <= 1e-9) continue;
    const sx = cx + (dx * rx + dy * ry + dz * rz) * focal / vz;
    if (sx < 0 || sx > W) continue;
    const sy = cy - (dx * ux + dy * uy + dz * uz) * focal / vz;
    if (sy < 0 || sy > H) continue;
    drawn++;
    const big = S.H[k] < 1000;                       /* H < 10: the large ones */
    const bucket = c * 2 + (big ? 1 : 0);
    if (bucket !== cur) { cur = bucket; ctx.fillStyle = 'rgba(' + SBDB_CLASS[c].rgb + ',' + ((big ? 0.85 : 0.35) * zoomFade[c]).toFixed(3) + ')'; }
    const s = big ? 2 : 1;
    ctx.fillRect(sx, sy, s, s);
    /* comets near the Sun grow a tail pointing away from it */
    if (c >= 4 && sunScreen) {
      const r = Math.hypot(S.px[k], S.py[k], S.pz[k]);
      if (r < 3.5) {
        let tx = sx - sunScreen.x, ty = sy - sunScreen.y;
        const m = Math.hypot(tx, ty) || 1; tx /= m; ty /= m;
        const len = clamp(22 / (r * r), 4, 60);
        const g = ctx.createLinearGradient(sx, sy, sx + tx * len, sy + ty * len);
        g.addColorStop(0, 'rgba(' + SBDB_CLASS[c].rgb + ',0.7)');
        g.addColorStop(1, 'rgba(' + SBDB_CLASS[c].rgb + ',0)');
        ctx.strokeStyle = g; ctx.lineWidth = 1.2;
        ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(sx + tx * len, sy + ty * len); ctx.stroke();
        cur = -1;
      }
    }
  }
  ctx.restore();
}

/* Promote a catalogue row into a real object. Elements are converted to the
   JPL-planet style keplerPos() expects (mean longitude at J2000, rate per
   century), so it rides the same orbit machinery as the planets. */
function bulkSbdbObject(k) {
  const id = 'sb' + k;
  if (byId[id]) return byId[id];
  const S = bulkSbdb, name = S.names[k], c = S.cls[k], cls = SBDB_CLASS[c];
  const nDeg = S.nrad[k] / DEG;                       /* deg per day */
  const ANG = 360 / 65536;
  const iDeg = S.I[k] * ANG, omDeg = S.OM[k] * ANG, wDeg = S.W[k] * ANG, mDeg = S.m0[k] / DEG;
  const varpi = wDeg + omDeg;
  /* mean longitude at J2000 = (M + varpi) at SBDB_EPOCH, wound back at the mean rate */
  const L0 = mDeg + varpi - nDeg * (S.epoch - J2000);
  const aAU = S.a[k], e = S.e[k];
  const isComet = c >= 4;
  const diamKm = S.diam[k] ? S.diam[k] / 10 : 1329 / Math.sqrt(0.15) * Math.pow(10, -S.H[k] / 500);
  const periodYr = Math.pow(aAU, 1.5);
  /* asteroids and comets are often named after the same person (2688 Halley vs 1P/Halley),
     so only take a featured description when the class agrees with the list it came from */
  const list = (typeof FEATURED_ASTEROIDS !== 'undefined') ? (isComet ? FEATURED_COMETS : FEATURED_ASTEROIDS) : [];
  const key = list.indexOf(S.keys[k]) >= 0 ? S.keys[k] : (list.indexOf(S.names[k]) >= 0 ? S.names[k] : null);
  const famous = key && FAMOUS_SMALL_BODIES[key];
  const o = addObject({
    id: id, name: name, type: cls.label, kind: 'dwarf', category: isComet ? 'comets' : 'asteroids',
    color: cls.color, color2: cls.color, radiusAU: (diamKm / 2) / AU_KM,
    el: { a: aAU, e: e, i: iDeg, L: L0, w: varpi, O: omDeg, aD: 0, eD: 0, iD: 0, LD: nDeg * 36525, wD: 0, OD: 0 },
    pole: { x: 0, y: 0, z: 1 }, prio: 58, minPx: 1.8,
    data: {
      desc: famous || ('From the JPL Small-Body Database. ' + cls.label + (S.diam[k] ? '' : '; size estimated from brightness') + '.'),
      stats: {
        'Class': cls.label,
        'Semi-major axis': aAU.toFixed(3) + ' AU',
        'Eccentricity': e.toFixed(3),
        'Inclination': iDeg.toFixed(1) + '°',
        'Perihelion': (aAU * (1 - e)).toFixed(3) + ' AU',
        'Orbital period': periodYr < 2 ? (periodYr * 365.25).toFixed(0) + ' days' : periodYr.toFixed(1) + ' years',
        'Diameter': S.diam[k] ? diamKm.toFixed(1) + ' km' : '≈' + diamKm.toFixed(1) + ' km (from H)',
        'Absolute magnitude H': (S.H[k] / 100).toFixed(1)
      }
    }
  });
  return o;
}


/* ---- Galaxies with distances (HyperLEDA). Layout, 32-bit first:
     f32 x,y,z (AU, ecliptic) | f32 radius (AU) | i16 mag*100 | i8 type (HyperLEDA t, -99 unknown)
     | u8 pa (0..255 = 0..180 deg) | u8 axis ratio *255.  Brightest first. */
let bulkGal = null;
if (typeof GAL_BIN !== 'undefined') {
  const b = b64buf(GAL_BIN), n = GAL_N;
  const X = new Float32Array(b, 0, n), Y = new Float32Array(b, 4 * n, n), Z = new Float32Array(b, 8 * n, n);
  const R = new Float32Array(b, 12 * n, n);
  const MG = new Int16Array(b, 16 * n, n), T = new Int8Array(b, 18 * n, n);
  const PA = new Uint8Array(b, 19 * n, n), AR = new Uint8Array(b, 20 * n, n);
  /* unit directions and distances, so the common case (camera inside the Milky
     Way, galaxies effectively at infinity) needs no per-frame subtraction */
  const dist = new Float32Array(n), dx = new Float32Array(n), dy = new Float32Array(n), dz = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const m = Math.hypot(X[i], Y[i], Z[i]) || 1; dist[i] = m; dx[i] = X[i] / m; dy[i] = Y[i] / m; dz[i] = Z[i] / m;
  }
  bulkGal = { n: n, x: X, y: Y, z: Z, r: R, mag: MG, t: T, pa: PA, ar: AR, dist: dist, dx: dx, dy: dy, dz: dz,
              names: GAL_NAMES.split('|') };
  bulkGal.lname = bulkGal.names.map(function (v) { return v.toLowerCase(); });
  { let lo = 0, hi = n; while (lo < hi) { const m = (lo + hi) >> 1; if (bulkGal.mag[m] <= 1350) lo = m + 1; else hi = m; } bulkGal.defaultN = lo; }
}
function galType(t) {
  if (t <= -90) return ['Galaxy', '#d0d6e8', 'galaxy'];
  if (t <= -3.5) return ['Elliptical galaxy', '#e6e0d2', 'galaxy'];
  if (t <= -0.5) return ['Lenticular galaxy', '#e2dccc', 'galaxy'];
  if (t < 8.5) return ['Spiral galaxy', '#cdd8f0', 'spiral'];
  return ['Irregular galaxy', '#b8d0ff', 'galaxy'];
}
/* screen direction of celestial north at a sky direction — used to orient the
   position angle of the few galaxies large enough to draw as ellipses */
const NORTH_ECL = eqToEcl(0, 0, 1);
function drawBulkGalaxies() {
  const G = bulkGal; if (!G || !opts.galaxies || (deep.on && deep.horizonFade <= 0.01)) return;
  const camR = Math.hypot(camPos.x, camPos.y, camPos.z);
  const far = camR > 1.5e6 * LY_AU;                /* well outside the Milky Way: real parallax */
  const px = camPos.x, py = camPos.y, pz = camPos.z;
  const fx = fwd.x, fy = fwd.y, fz = fwd.z, rx = right.x, ry = right.y, rz = right.z, ux = up.x, uy = up.y, uz = up.z;
  let drawn = 0, grads = 0, bucket = -1;
  const budget = 40000, GRAD = 320;
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  const nDraw = opts.allGalaxies ? G.n : G.defaultN;
  const farFade = (1 - 0.85 * clamp((camR / LY_AU - 3e9) / 2e10, 0, 1)) * (deep.on ? deep.horizonFade : 1);   /* a speck from the edge of the universe */
  for (let i = 0; i < nDraw && drawn < budget; i++) {
    let vz, sx, sy, range;
    if (far) {
      const ddx = G.x[i] - px, ddy = G.y[i] - py, ddz = G.z[i] - pz;
      vz = ddx * fx + ddy * fy + ddz * fz; if (vz <= 1e-6) continue;
      sx = cx + (ddx * rx + ddy * ry + ddz * rz) * focal / vz; if (sx < -30 || sx > W + 30) continue;
      sy = cy - (ddx * ux + ddy * uy + ddz * uz) * focal / vz; if (sy < -30 || sy > H + 30) continue;
      range = Math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
    } else {
      const ddx = G.dx[i], ddy = G.dy[i], ddz = G.dz[i];
      vz = ddx * fx + ddy * fy + ddz * fz; if (vz <= 0.02) continue;
      sx = cx + (ddx * rx + ddy * ry + ddz * rz) * focal / vz; if (sx < -30 || sx > W + 30) continue;
      sy = cy - (ddx * ux + ddy * uy + ddz * uz) * focal / vz; if (sy < -30 || sy > H + 30) continue;
      range = G.dist[i];
    }
    drawn++;
    const m = G.mag[i] * 0.01;
    const rPx = (G.r[i] / range) * focal;
    const gt = galType(G.t[i]);
    if (rPx < 2.4 || grads >= GRAD) {
      const b = (m < 11 ? 0 : m < 13 ? 1 : m < 14.5 ? 2 : 3) + (gt[2] === 'spiral' ? 0 : 4);
      if (b !== bucket) { bucket = b; ctx.fillStyle = hexA(gt[1], [0.55, 0.36, 0.20, 0.10][b & 3] * farFade); }
      const s = rPx < 1.2 ? 1 : Math.min(3, rPx);
      ctx.fillRect(sx - s * 0.5, sy - s * 0.5, s, s);
    } else {
      grads++; bucket = -1;
      /* orientation: position angle is measured from north through east on the sky */
      const nq = projectDir({ x: G.dx[i] + NORTH_ECL.x * 1e-3, y: G.dy[i] + NORTH_ECL.y * 1e-3, z: G.dz[i] + NORTH_ECL.z * 1e-3 });
      const ang = nq ? Math.atan2(nq.y - sy, nq.x - sx) : -Math.PI / 2;
      const pa = G.pa[i] / 255 * Math.PI;
      ctx.save();
      ctx.translate(sx, sy); ctx.rotate(ang + pa); ctx.scale(1, Math.max(0.15, G.ar[i] / 255));
      const g = ctx.createRadialGradient(0, 0, 0, 0, 0, rPx);
      g.addColorStop(0, 'rgba(255,255,255,0.55)');
      g.addColorStop(0.15, hexA(gt[1], 0.36));
      g.addColorStop(0.5, hexA(gt[1], 0.12));
      g.addColorStop(1, hexA(gt[1], 0));
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(0, 0, rPx, 0, TAU); ctx.fill();
      ctx.restore();
    }
  }
  ctx.restore();
}
function bulkGalObject(i) {
  const id = 'gal' + i;
  if (byId[id]) return byId[id];
  const G = bulkGal, gt = galType(G.t[i]);
  const ly = G.dist[i] / LY_AU, name = G.names[i];
  const o = addObject({
    id: id, name: name, type: gt[0], kind: 'deepsky', sub: gt[2], category: 'galaxies',
    color: gt[1], color2: gt[1], radiusAU: G.r[i],
    fixed: { x: G.x[i], y: G.y[i], z: G.z[i] }, prio: 40, minPx: 2.0,
    data: { d: ly, tilt: G.pa[i] / 255 * Math.PI, flat: Math.max(0.15, G.ar[i] / 255),
            desc: 'From the HyperLEDA extragalactic database. ' + gt[0] + (ly > 1e6 ? ', about ' + (ly / 1e6).toFixed(1) + ' million light years away.' : '.'),
            stats: { 'Type': gt[0], 'Distance': fmtLy(ly), 'Diameter': fmtLy(G.r[i] * 2 / LY_AU),
                     'Magnitude (B)': (G.mag[i] * 0.01).toFixed(1), 'Catalogue': name.indexOf('PGC') === 0 ? name : 'PGC via HyperLEDA' } }
  });
  o.pos = o.fixed;
  return o;
}

/* ---- Artificial satellites (Celestrak GP elements), Earth-centred. Layout, 32-bit first:
     f32 a (km) | u16 e*65535 | u16 i | u16 raan | u16 argp | u16 M0 (0..65535 = 0..360)
     | u8 class.  Propagated to the common epoch SAT_EPOCH at build time. */
let bulkSat = null;
if (typeof SAT_BIN !== 'undefined') {
  const b = b64buf(SAT_BIN), n = SAT_N;
  const A = new Float32Array(b, 0, n);
  const E = new Uint16Array(b, 4 * n, n), I = new Uint16Array(b, 6 * n, n), RA = new Uint16Array(b, 8 * n, n);
  const AP = new Uint16Array(b, 10 * n, n), M0 = new Uint16Array(b, 12 * n, n), CL = new Uint8Array(b, 14 * n, n);
  const ANG = TAU / 65536, MU = 398600.4418;         /* km^3 s^-2 */
  const nrad = new Float32Array(n), ecc = new Float32Array(n), m0 = new Float32Array(n);
  const Px = new Float32Array(n), Py = new Float32Array(n), Pz = new Float32Array(n);
  const Qx = new Float32Array(n), Qy = new Float32Array(n), Qz = new Float32Array(n);
  for (let k = 0; k < n; k++) {
    const e = E[k] / 65535, inc = I[k] * ANG, O = RA[k] * ANG, w = AP[k] * ANG;
    const cw = Math.cos(w), sw = Math.sin(w), cO = Math.cos(O), sO = Math.sin(O), cI = Math.cos(inc), sI = Math.sin(inc);
    Px[k] = cw * cO - sw * sO * cI;  Py[k] = cw * sO + sw * cO * cI;  Pz[k] = sw * sI;
    Qx[k] = -sw * cO - cw * sO * cI; Qy[k] = -sw * sO + cw * cO * cI; Qz[k] = cw * sI;
    ecc[k] = e; m0[k] = M0[k] * ANG; nrad[k] = Math.sqrt(MU / (A[k] * A[k] * A[k])) * 86400;   /* rad/day */
  }
  bulkSat = { n: n, a: A, e: ecc, m0: m0, nrad: nrad, Px, Py, Pz, Qx, Qy, Qz, cls: CL, I: I, RA: RA, AP: AP,
              px: new Float32Array(n), py: new Float32Array(n), pz: new Float32Array(n),
              names: SAT_NAMES.split('|'), epoch: SAT_EPOCH, lastJd: -1 };
  bulkSat.lname = bulkSat.names.map(function (v) { return v.toLowerCase(); });
}
const SAT_CLASS = [
  { label: 'Space station',           color: '#ffffff', rgb: '255,255,255' },
  { label: 'Starlink satellite',      color: '#8fd3ff', rgb: '143,211,255' },
  { label: 'OneWeb satellite',        color: '#b8a9ff', rgb: '184,169,255' },
  { label: 'Navigation satellite',    color: '#ffe08a', rgb: '255,224,138' },
  { label: 'Geostationary satellite', color: '#ffb072', rgb: '255,176,114' },
  { label: 'Satellite',               color: '#c8d0dc', rgb: '200,208,220' }
];
/* ECI (equatorial, km) -> ecliptic AU, added to Earth's position */
const SAT_CE = Math.cos(OBLIQ), SAT_SE = Math.sin(OBLIQ);
function satSolveOne(S, k, jd) {
  const e0 = byId.earth.pos, CE = SAT_CE, SE = SAT_SE;
  {
    const e = S.e[k];
    let M = S.m0[k] + S.nrad[k] * (jd - S.epoch); M -= TAU * Math.floor(M / TAU);
    let Ea = M;
    for (let it = 0; it < 6; it++) { const d = (Ea - e * Math.sin(Ea) - M) / (1 - e * Math.cos(Ea)); Ea -= d; if (Math.abs(d) < 1e-6) break; }
    const a = S.a[k], xp = a * (Math.cos(Ea) - e), yp = a * Math.sqrt(1 - e * e) * Math.sin(Ea);
    const ex = (xp * S.Px[k] + yp * S.Qx[k]) / AU_KM, ey = (xp * S.Py[k] + yp * S.Qy[k]) / AU_KM, ez = (xp * S.Pz[k] + yp * S.Qz[k]) / AU_KM;
    S.px[k] = e0.x + ex;
    S.py[k] = e0.y + ey * CE + ez * SE;
    S.pz[k] = e0.z - ey * SE + ez * CE;
  }
}
/* Solving 16,500 orbits every frame cost ~12 ms. A third per frame is invisible at real
   time (a satellite moves ~1 deg between refreshes) and fast time is a blur regardless. */
function satUpdate(jd) {
  const S = bulkSat; if (!S) return;
  if (S.lastJd < 0 || Math.abs(jd - S.lastJd) > 0.5) {
    for (let k = 0; k < S.n; k++) satSolveOne(S, k, jd);
    S.lastJd = jd; S.cursor = 0; return;
  }
  const step = Math.ceil(S.n / 3), end = Math.min(S.n, (S.cursor || 0) + step);
  for (let k = S.cursor || 0; k < end; k++) satSolveOne(S, k, jd);
  S.cursor = end >= S.n ? 0 : end;
  if (S.cursor === 0) S.lastJd = jd;
}
function satsVisible() {
  if (!bulkSat || !opts.satellites) return false;
  /* only while Earth or one of its satellites is the subject — never as ambient clutter */
  const sel = selectedId && byId[selectedId];
  const inContext = (sel && (sel.id === 'earth' || sel.category === 'satellites')) || focusCategory === 'satellites';
  if (!inContext) return false;
  const e = byId.earth.pos;
  return Math.hypot(camPos.x - e.x, camPos.y - e.y, camPos.z - e.z) < 0.02;   /* within ~3 million km */
}
function drawSatellites() {
  if (!satsVisible() || deep.on) return;
  const S = bulkSat;
  satUpdate(simJd);
  const px = camPos.x, py = camPos.y, pz = camPos.z;
  const fx = fwd.x, fy = fwd.y, fz = fwd.z, rx = right.x, ry = right.y, rz = right.z, ux = up.x, uy = up.y, uz = up.z;
  let cur = -1;
  for (let k = 0; k < S.n; k++) {
    const dx = S.px[k] - px, dy = S.py[k] - py, dz = S.pz[k] - pz;
    const vz = dx * fx + dy * fy + dz * fz; if (vz <= 1e-12) continue;
    const sx = cx + (dx * rx + dy * ry + dz * rz) * focal / vz; if (sx < 0 || sx > W) continue;
    const sy = cy - (dx * ux + dy * uy + dz * uz) * focal / vz; if (sy < 0 || sy > H) continue;
    const c = S.cls[k];
    if (c !== cur) { cur = c; ctx.fillStyle = 'rgba(' + SAT_CLASS[c].rgb + ',' + (c === 0 ? 1 : (c === 3 || c === 4) ? 0.85 : 0.42) + ')'; }
    const s = c === 0 ? 3 : 1.4;
    ctx.fillRect(sx - s * 0.5, sy - s * 0.5, s, s);
  }
}
function bulkSatObject(k) {
  const id = 'sat' + k;
  if (byId[id]) return byId[id];
  const S = bulkSat, cls = SAT_CLASS[S.cls[k]], name = S.names[k];
  const ANG = 360 / 65536;
  const a = S.a[k], e = S.e[k], perMin = TAU / S.nrad[k] * 1440;
  let famous = null;
  if (typeof FAMOUS_SATELLITES !== 'undefined') {
    let best = '';
    for (const key in FAMOUS_SATELLITES) if (name.toUpperCase().indexOf(key.toUpperCase()) === 0 && key.length > best.length) best = key;
    if (best) famous = FAMOUS_SATELLITES[best];
  }
  const o = addObject({
    id: id, name: name, type: cls.label, kind: 'probe', category: 'satellites', color: cls.color,
    radiusAU: 0, mode: 'gp', satIndex: k, prio: 63, minPx: 0,
    data: { desc: famous || ('From Celestrak\'s active-satellite catalogue. Position is a two-body propagation of the published mean elements, accurate near the element epoch.'),
            stats: { 'Class': cls.label, 'Altitude (perigee)': Math.round(a * (1 - e) - 6371).toLocaleString() + ' km',
                     'Altitude (apogee)': Math.round(a * (1 + e) - 6371).toLocaleString() + ' km',
                     'Orbital period': perMin < 120 ? perMin.toFixed(1) + ' min' : (perMin / 60).toFixed(2) + ' hours',
                     'Inclination': (S.I[k] * ANG).toFixed(1) + '°', 'Eccentricity': e.toFixed(4) } }
  });
  return o;
}

/* colour from B-V index, in coarse steps */
function ciColor(ci) {
  if (ci < -0.10) return '#a8c0ff';
  if (ci < 0.15) return '#cad8ff';
  if (ci < 0.45) return '#f2f4ff';
  if (ci < 0.80) return '#fff6e8';
  if (ci < 1.30) return '#ffd9a8';
  return '#ffb27a';
}

/* ============================================================
   4. Positions for a given time
   ============================================================ */
let simJd = jdFromDate(new Date());

/* Position of an L2 observatory: on the Sun-Earth line 1.5 million km
   beyond Earth, reached by a ~30 day cruise, then held in a halo orbit
   around the point rather than sitting exactly on it. */
function l2Position(o, jd) {
  const e = keplerPos(byId.earth.el, (jd - J2000) / 36525);
  const m = Math.hypot(e.x, e.y, e.z) || 1;
  const dir = { x: e.x / m, y: e.y / m, z: e.z / m };
  let u = cross({ x: 0, y: 0, z: 1 }, dir);
  if (Math.hypot(u.x, u.y, u.z) < 1e-9) u = { x: 1, y: 0, z: 0 };
  u = norm(u);
  const v = cross(dir, u);

  let f = 1, cruising = false;
  if (o.launchJd != null) {
    if (jd <= o.launchJd) { f = 0; cruising = true; }
    else if (jd < o.arriveJd) { f = (jd - o.launchJd) / (o.arriveJd - o.launchJd); cruising = true; }
  }
  const k = (m + o.offsetAU * f) / m;
  const p = { x: e.x * k, y: e.y * k, z: e.z * k, cruising: cruising, f: f };
  if (cruising) {
    /* a transfer curves away from the straight line and settles back */
    const lat = Math.sin(Math.PI * f) * (o.haloAU || 0.004) * 0.9;
    p.x += u.x * lat; p.y += u.y * lat; p.z += u.z * lat;
  } else if (o.haloAU) {
    const th = TAU * ((jd - J2000) / (o.haloDays || 180));
    const c = Math.cos(th) * o.haloAU, sn = Math.sin(th) * o.haloAU * 0.55;
    p.x += u.x * c + v.x * sn;
    p.y += u.y * c + v.y * sn;
    p.z += u.z * c + v.z * sn;
  }
  return p;
}

function updatePositions(jd) {
  const T = (jd - J2000) / 36525;

  for (let i = 0; i < objects.length; i++) {
    const o = objects[i];
    if (o.fixed) { o.pos = o.fixed; continue; }
    if (o.kind === 'sun') { o.pos.x = o.pos.y = o.pos.z = 0; continue; }
    if (o.el) {
      const p = keplerPos(o.el, T);
      o.pos.x = p.x; o.pos.y = p.y; o.pos.z = p.z;
    }
  }
  sbdbUpdate(jd);

  /* moons, after their parents */
  for (let i = 0; i < objects.length; i++) {
    const o = objects[i];
    if (o.kind !== 'moon') continue;
    const par = byId[o.parent];
    if (o.id === 'moon') {
      /* Our Moon gets a real ephemeris (Meeus' low-precision series, ~0.3°), so its
         phase, its place in the sky and the timing of full and new Moon are genuine. */
      const T = (jd - J2000) / 36525;
      const Lp = (218.3165 + 481267.8813 * T) * DEG, Mp = (134.963 + 477198.8676 * T) * DEG;
      const D = (297.850 + 445267.1115 * T) * DEG, M = (357.529 + 35999.0503 * T) * DEG, F = (93.272 + 483202.0175 * T) * DEG;
      const lon = Lp + (6.289 * Math.sin(Mp) + 1.274 * Math.sin(2 * D - Mp) + 0.658 * Math.sin(2 * D) + 0.214 * Math.sin(2 * Mp) - 0.186 * Math.sin(M) - 0.114 * Math.sin(2 * F)) * DEG;
      const lat = (5.128 * Math.sin(F) + 0.281 * Math.sin(Mp + F) - 0.278 * Math.sin(F - Mp) - 0.173 * Math.sin(F - 2 * D)) * DEG;
      const dist = (385001 - 20905 * Math.cos(Mp) - 3699 * Math.cos(2 * D - Mp) - 2956 * Math.cos(2 * D)) / AU_KM;
      const cl = Math.cos(lat);
      o.pos.x = par.pos.x + dist * cl * Math.cos(lon);
      o.pos.y = par.pos.y + dist * cl * Math.sin(lon);
      o.pos.z = par.pos.z + dist * Math.sin(lat);
      continue;
    }
    /* a moon's inclination is quoted against its parent's equator, except where the
       catalogue says otherwise (the Moon's 5.1 deg is from the ecliptic) */
    const p = o.refPlane === 'ecliptic' ? { x: 0, y: 0, z: 1 } : par.pole;
    let u = cross({ x: 0, y: 0, z: 1 }, p);
    if (Math.hypot(u.x, u.y, u.z) < 1e-6) u = { x: 1, y: 0, z: 0 };
    u = norm(u);
    const v = cross(p, u);
    const th = TAU * ((jd - J2000) / o.periodDays);
    const ci = Math.cos(o.inc), si = Math.sin(o.inc);
    const c = Math.cos(th) * o.aAU, s = Math.sin(th) * o.aAU;
    o.pos.x = par.pos.x + c * u.x + s * (v.x * ci + p.x * si);
    o.pos.y = par.pos.y + c * u.y + s * (v.y * ci + p.y * si);
    o.pos.z = par.pos.z + c * u.z + s * (v.z * ci + p.z * si);
  }
  /* exoplanets, around their host stars */
  for (let i = 0; i < objects.length; i++) {
    const o = objects[i];
    /* bulk catalogue rows are placed at a fixed system position and carry no host
       object — dereferencing one here threw and permanently killed the frame loop */
    if (o.kind !== 'exoplanet' || !o.parent) continue;
    const host = byId[o.parent];
    if (!host) continue;
    const th = TAU * ((jd - J2000) / o.periodDays);
    const c = Math.cos(th) * o.aAU, sn = Math.sin(th) * o.aAU;
    o.pos.x = host.pos.x + c * o.ou.x + sn * o.ov.x;
    o.pos.y = host.pos.y + c * o.ou.y + sn * o.ov.y;
    o.pos.z = host.pos.z + c * o.ou.z + sn * o.ov.z;
  }

  /* spacecraft */
  for (let i = 0; i < objects.length; i++) {
    const o = objects[i];
    if (o.kind !== 'probe') continue;
    if (o.mode === 'escape') {
      const d = Math.max(1, o.distAU + o.speedAUday * (jd - o.epochJd));
      o.pos.x = o.dir.x * d; o.pos.y = o.dir.y * d; o.pos.z = o.dir.z * d;
      o.liveDist = d;
    } else if (o.mode === 'gp') {
      const S = bulkSat;
      if (S) { satSolveOne(S, o.satIndex, jd); o.pos.x = S.px[o.satIndex]; o.pos.y = S.py[o.satIndex]; o.pos.z = S.pz[o.satIndex]; }
    } else if (o.mode === 'earthorbit' || o.mode === 'orbiter') {
      /* a satellite of another body — same geometry as a moon, with an optional
         eccentric orbit for the long thin science orbits (Juno, Akatsuki) */
      const par = byId[o.parent] || byId.earth, pl = par.pole || { x: 0, y: 0, z: 1 };
      let u = cross({ x: 0, y: 0, z: 1 }, pl);
      if (Math.hypot(u.x, u.y, u.z) < 1e-6) u = { x: 1, y: 0, z: 0 };
      u = norm(u);
      const v = cross(pl, u);
      const e = o.ecc || 0;
      let M = TAU * ((jd - J2000) / o.periodDays); M -= TAU * Math.floor(M / TAU);
      let Ea = M;
      for (let it = 0; it < 8; it++) { const dd = (Ea - e * Math.sin(Ea) - M) / (1 - e * Math.cos(Ea)); Ea -= dd; if (Math.abs(dd) < 1e-8) break; }
      const xp = o.aAU * (Math.cos(Ea) - e), yp = o.aAU * Math.sqrt(1 - e * e) * Math.sin(Ea);
      const cw = Math.cos(o.argPeri || 0), sw = Math.sin(o.argPeri || 0);
      const c = xp * cw - yp * sw, sn = xp * sw + yp * cw;
      const ci = Math.cos(o.inc), si = Math.sin(o.inc);
      o.pos.x = par.pos.x + c * u.x + sn * (v.x * ci + pl.x * si);
      o.pos.y = par.pos.y + c * u.y + sn * (v.y * ci + pl.y * si);
      o.pos.z = par.pos.z + c * u.z + sn * (v.z * ci + pl.z * si);
    } else if (o.mode === 'atbody') {
      /* a lander, or a mission that ended in the body: a marker just off the sunlit limb */
      const par = byId[o.parent];
      const m = Math.hypot(par.pos.x, par.pos.y, par.pos.z) || 1, k = par.radiusAU * 1.25;
      o.pos.x = par.pos.x - par.pos.x / m * k;
      o.pos.y = par.pos.y - par.pos.y / m * k;
      o.pos.z = par.pos.z - par.pos.z / m * k;
    } else if (o.mode === 'trailing') {
      /* Earth-trailing: same orbit as Earth, falling further behind each year */
      const yrs = Math.max(0, (jd - o.launchJd) / 365.25);
      const lagDays = (o.driftDegYr * yrs / 360) * 365.25;
      const q = keplerPos(byId.earth.el, (jd - lagDays - J2000) / 36525);
      o.pos.x = q.x; o.pos.y = q.y; o.pos.z = q.z;
      o.lagDeg = o.driftDegYr * yrs;
    } else if (o.mode === 'l2' || o.mode === 'l2halo') {
      const q = l2Position(o, jd);
      o.pos.x = q.x; o.pos.y = q.y; o.pos.z = q.z;
      o.cruising = q.cruising;
      o.cruiseFrac = q.f;
    }
  }
}

/* ============================================================
   5. Camera
   ============================================================ */
const cam = {
  tx: 0, ty: 0, tz: 0,              // current target
  gx: 0, gy: 0, gz: 0,              // goal target
  dist: 48, distGoal: 48,
  yaw: -1.05, yawGoal: -1.05,
  pitch: 0.42, pitchGoal: 0.42,
  roll: 0, rollGoal: 0,             // tilt of the horizon; stays 0 unless free rotation is on
  pendYaw: 0, pendPitch: 0,         // free rotation still to be applied (smoothed like the goals)
  fov: 50 * DEG, fovGoal: 50 * DEG,
  follow: null,
  flight: null
};
let camPos = { x: 0, y: 0, z: 0 }, fwd = { x: 0, y: 0, z: 1 }, right = { x: 1, y: 0, z: 0 }, up = { x: 0, y: 1, z: 0 };
let focal = 600, cx = 0, cy = 0;

/* Deep time: the future of the universe, played forward on a log clock (see the
   module near the end of this file). Declared here so the draw passes can read it. */
const deep = { on: false, logY: 0, Y: 1, playing: false, hold: 0, chapter: -1, takeover: false,
               yrsPerSec: 0, horizonFade: 1, galaxyFade: 1, merge: 0, sunMaxR: 1, clock: 0 };
const MIN_DIST = 1e-5, MAX_DIST = 1e16;        /* out to beyond the edge of the observable universe */
/* The sky can be magnified like a telescope: the field of view narrows from the
   50° walk-around view down to about an arcminute, enough to fill the screen with
   a deep field that is a fraction of a pixel wide at normal zoom. */
const FOV_DEFAULT = 50 * DEG, FOV_MIN = 1.2 / 60 * DEG;

function updateCamera(dt) {
  /* fly-to animation */
  if (cam.flight) {
    const f = cam.flight;
    f.t += dt / f.dur;
    const t = clamp(f.t, 0, 1);
    const e = t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;   // easeInOutCubic
    const dest = f.obj;
    cam.gx = lerp(f.fx, dest.pos.x, e);
    cam.gy = lerp(f.fy, dest.pos.y, e);
    cam.gz = lerp(f.fz, dest.pos.z, e);
    cam.distGoal = Math.exp(lerp(Math.log(f.fd), Math.log(f.td), e));
    if (f.yaw != null) cam.yawGoal = lerp(f.fyaw, f.yaw, e);
    if (f.pitch != null) cam.pitchGoal = lerp(f.fpitch, f.pitch, e);
    cam.rollGoal = 0; cam.pendYaw = cam.pendPitch = 0;       /* a flight always arrives level */
    if (t >= 1) { cam.follow = f.noFollow ? null : dest.id; cam.flight = null; updateFocusChip(); }
  }
  const following = cam.follow && byId[cam.follow] && !cam.flight;
  if (following) {
    const o = byId[cam.follow];
    cam.gx = o.pos.x; cam.gy = o.pos.y; cam.gz = o.pos.z;
  }

  const k = 1 - Math.pow(0.0016, dt);       // frame-rate independent damping
  if (following) {
    /* Lock rigidly onto the body — no damping. A damped target always lags,
       and up close that is fatal: with the clock at 1 day/s Pluto covers
       ~400,000 km per real second while the camera orbits 8,000 km out, so
       the body outruns the target and leaves you staring at empty space. */
    cam.tx = cam.gx; cam.ty = cam.gy; cam.tz = cam.gz;
  } else {
    cam.tx = lerp(cam.tx, cam.gx, k);
    cam.ty = lerp(cam.ty, cam.gy, k);
    cam.tz = lerp(cam.tz, cam.gz, k);
  }
  if (!opts.freeRotate && !deep.on) { cam.rollGoal = 0; cam.pendYaw = cam.pendPitch = 0; }   /* the level horizon returns */
  cam.yaw = lerp(cam.yaw, cam.yawGoal, k);
  cam.pitch = lerp(cam.pitch, cam.pitchGoal, k);
  cam.roll = lerp(cam.roll, cam.rollGoal, k);
  cam.dist = Math.exp(lerp(Math.log(cam.dist), Math.log(cam.distGoal), k));
  cam.fov = Math.exp(lerp(Math.log(cam.fov), Math.log(cam.fovGoal), k));
  focal = (H / 2) / Math.tan(cam.fov / 2);

  /* Free rotation: a trackball. The drag turns the view about the camera's own
     axes, so there is no pole to stop at; the horizon is free to tilt. */
  if (opts.freeRotate && !cam.flight && Math.abs(cam.pendYaw) + Math.abs(cam.pendPitch) > 1e-7) {
    const a = cam.pendYaw * k, b = cam.pendPitch * k;
    cam.pendYaw -= a; cam.pendPitch -= b;
    const F0 = camFrame(cam.yaw, cam.pitch, cam.roll);
    const d1 = rotateAbout(F0.d, F0.up, a);                       /* sideways: about the camera's up */
    const cb = Math.cos(b), sb = Math.sin(b);
    const d2 = { x: d1.x * cb + F0.up.x * sb, y: d1.y * cb + F0.up.y * sb, z: d1.z * cb + F0.up.z * sb };
    const u2 = { x: F0.up.x * cb - d1.x * sb, y: F0.up.y * cb - d1.y * sb, z: F0.up.z * cb - d1.z * sb };
    const e = camAngles(d2, u2, cam.yaw);
    cam.yaw = cam.yawGoal = e.yaw; cam.pitch = cam.pitchGoal = e.pitch; cam.roll = cam.rollGoal = e.roll;
  }
  const Fr = camFrame(cam.yaw, cam.pitch, cam.roll);
  if (cam.fov < FOV_DEFAULT * 0.98 && byId.earth) {
    /* Telescope view: stand on Earth and look outward. The orbit target is pushed
       out along the line of sight so the camera position lands exactly on Earth,
       instead of the Sun and planets sitting in front of the patch of sky. */
    const e = byId.earth.pos;
    cam.tx = e.x - cam.dist * Fr.d.x;
    cam.ty = e.y - cam.dist * Fr.d.y;
    cam.tz = e.z - cam.dist * Fr.d.z;
  }
  camPos = { x: cam.tx + cam.dist * Fr.d.x, y: cam.ty + cam.dist * Fr.d.y, z: cam.tz + cam.dist * Fr.d.z };
  fwd = { x: -Fr.d.x, y: -Fr.d.y, z: -Fr.d.z };
  right = Fr.right;
  up = Fr.up;
}
/* The camera's axes from its three angles. With no roll this is the classic
   turntable: north stays up. The closed form has no special case at the poles. */
function camFrame(yaw, pitch, roll) {
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const d = { x: cp * cy, y: cp * sy, z: sp };                  /* from the target to the camera */
  const r0 = { x: -sy, y: cy, z: 0 }, u0 = { x: -sp * cy, y: -sp * sy, z: cp };
  if (!roll) return { d: d, right: r0, up: u0 };
  const cr = Math.cos(roll), sr = Math.sin(roll);
  return { d: d,
           right: { x: r0.x * cr + u0.x * sr, y: r0.y * cr + u0.y * sr, z: r0.z * cr + u0.z * sr },
           up: { x: u0.x * cr - r0.x * sr, y: u0.y * cr - r0.y * sr, z: u0.z * cr - r0.z * sr } };
}
/* ...and back: the three angles that describe a given view direction and up vector */
function camAngles(d, upv, prevYaw) {
  const pitch = Math.asin(clamp(d.z, -1, 1));
  let yaw = Math.hypot(d.x, d.y) > 1e-6 ? Math.atan2(d.y, d.x) : prevYaw;
  yaw += TAU * Math.round((prevYaw - yaw) / TAU);
  const L = camFrame(yaw, pitch, 0);
  const roll = Math.atan2(-(upv.x * L.right.x + upv.y * L.right.y + upv.z * L.right.z), upv.x * L.up.x + upv.y * L.up.y + upv.z * L.up.z);
  return { yaw: yaw, pitch: pitch, roll: roll };
}
function rotateAbout(v, axis, ang) {
  const c = Math.cos(ang), s = Math.sin(ang), dot = v.x * axis.x + v.y * axis.y + v.z * axis.z, x = cross(axis, v);
  return { x: v.x * c + x.x * s + axis.x * dot * (1 - c), y: v.y * c + x.y * s + axis.y * dot * (1 - c), z: v.z * c + x.z * s + axis.z * dot * (1 - c) };
}

function project(p) {
  const dx = p.x - camPos.x, dy = p.y - camPos.y, dz = p.z - camPos.z;
  const vz = dx * fwd.x + dy * fwd.y + dz * fwd.z;
  if (vz <= 1e-11) return null;
  const vx = dx * right.x + dy * right.y + dz * right.z;
  const vy = dx * up.x + dy * up.y + dz * up.z;
  return { x: cx + vx * focal / vz, y: cy - vy * focal / vz, z: vz };
}
function projectDir(d) {
  const vz = d.x * fwd.x + d.y * fwd.y + d.z * fwd.z;
  if (vz <= 0.001) return null;
  const vx = d.x * right.x + d.y * right.y + d.z * right.z;
  const vy = d.x * up.x + d.y * up.y + d.z * up.z;
  return { x: cx + vx * focal / vz, y: cy - vy * focal / vz, z: vz };
}

/* ============================================================
   6. Renderer
   ============================================================ */
const canvas = document.getElementById('scene');
const ctx = canvas.getContext('2d', { alpha: false });
const stage = document.getElementById('stage');
let W = 0, H = 0, DPR = 1;

function resize() {
  /* Ignore zero-size layouts (hidden panes, pre-layout first paint) — resizing
     to 0 would put the projection centre at 0,0 and throw everything into the
     top-left corner until the next real resize. */
  if (!stage.clientWidth || !stage.clientHeight) return;
  DPR = Math.min(window.devicePixelRatio || 1, 2);
  W = stage.clientWidth; H = stage.clientHeight;
  canvas.width = Math.round(W * DPR);
  canvas.height = Math.round(H * DPR);
  canvas.style.width = W + 'px';
  canvas.style.height = H + 'px';
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  cx = W / 2; cy = H / 2;
  focal = (H / 2) / Math.tan(cam.fov / 2);
}
window.addEventListener('resize', resize);
/* The stage can change size without a window resize (panes, split views, a
   layout that settles after first paint). Watch the element itself, so the
   projection centre is never left at 0,0. */
if (window.ResizeObserver) {
  const ro = new ResizeObserver(function () {
    if (stage.clientWidth !== W || stage.clientHeight !== H) resize();
  });
  ro.observe(stage);
}

const opts = {
  orbits: true, labels: true, moons: true, belts: true, oort: true,
  stars: true, ngc: true, milkyway: true,
  mainbelt: true, neo: true, trojans: true, tno: true, comets: true,
  galaxies: true, satellites: false,
  allStars: false, allGalaxies: false, allSmallBodies: false, probes: true, realSize: false,
  constellations: true, deepfields: true, structures: 'shells', surfaces: false, freeRotate: false,
  mwGain: 0.30           /* Milky Way brightness — calibrated against the reference */
};

let selectedId = null;
let focusCategory = null;   /* which category panel is open */
let labelBoxes = [], backdropLabelBoxes = [];

/* the painted sky is only right from inside our own neighbourhood */
function skyFade() {
  const camR = Math.hypot(camPos.x, camPos.y, camPos.z);
  return clamp(1 - (camR / LY_AU - 300) / 2000, 0, 1);
}
function overlapsLabel(box) {
  /* backdrop labels (regions, walls, the galaxy) live in their own list until the object pass copies them */
  const lists = [labelBoxes, backdropLabelBoxes];
  for (let l = 0; l < 2; l++) {
    const L = lists[l];
    for (let b = 0; b < L.length; b++) {
      const q = L[b];
      if (box.x < q.x + q.w && box.x + box.w > q.x && box.y < q.y + q.h && box.y + box.h > q.y) return true;
    }
  }
  return false;
}


/* The galaxy from outside. Fades in as the camera leaves the solar neighbourhood
   (the real catalogue stars and the painted band represent the inside view), and
   turns into a soft oriented ellipse — like every other galaxy in the map — when
   it is small on screen. */
function drawMilkyWay() {
  if (!MW || !opts.milkyway) return;
  const camR = Math.hypot(camPos.x, camPos.y, camPos.z) / LY_AU;
  const fade = clamp((camR - 1500) / 4500, 0, 1) * (deep.on ? deep.galaxyFade * (1 - 0.88 * deep.merge) : 1);
  if (fade <= 0.005) return;
  const C = MW.centre;
  const c = project(C);
  const distC = Math.hypot(C.x - camPos.x, C.y - camPos.y, C.z - camPos.z);
  const galPx = (MW.radius / distC) * focal;               /* disc radius on screen, roughly */
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';

  if (c && galPx < 2) {                                    /* a dot, like the catalogue galaxies at this range */
    ctx.fillStyle = 'rgba(230,225,245,' + (0.8 * fade).toFixed(3) + ')';
    ctx.fillRect(c.x - 1, c.y - 1, 2, 2);
    ctx.restore(); return;
  }
  if (c && galPx < 140) {                                  /* far: one oriented soft ellipse */
    const R = MW.radius;
    const pu = project({ x: C.x + MW.ex.x * R, y: C.y + MW.ex.y * R, z: C.z + MW.ex.z * R });
    const pv = project({ x: C.x + MW.ey.x * R, y: C.y + MW.ey.y * R, z: C.z + MW.ey.z * R });
    if (pu && pv) {
      const sa = clamp((140 - galPx) / 60, 0, 1) * fade;
      ctx.save();
      ctx.transform(pu.x - c.x, pu.y - c.y, pv.x - c.x, pv.y - c.y, c.x, c.y);
      const g = ctx.createRadialGradient(0, 0, 0, 0, 0, 1);
      g.addColorStop(0, 'rgba(255,252,245,' + (0.55 * sa).toFixed(3) + ')');
      g.addColorStop(0.14, 'rgba(240,228,205,' + (0.34 * sa).toFixed(3) + ')');
      g.addColorStop(0.5, 'rgba(195,205,242,' + (0.14 * sa).toFixed(3) + ')');
      g.addColorStop(1, 'rgba(180,190,240,0)');
      ctx.fillStyle = g;
      ctx.beginPath(); ctx.arc(0, 0, 1, 0, TAU); ctx.fill();
      ctx.restore();
    }
  }
  const detail = clamp((galPx - 40) / 60, 0, 1) * fade;
  if (detail > 0.01) {
    const fx = fwd.x, fy = fwd.y, fz = fwd.z, rx = right.x, ry = right.y, rz = right.z, ux = up.x, uy = up.y, uz = up.z;
    const px0 = camPos.x, py0 = camPos.y, pz0 = camPos.z;
    /* soft light and dust, far to near so dust in front darkens what is behind it */
    const B = MW.blobs, ord = MW.order, d2 = MW.dist2;
    for (let i = 0; i < B.length; i++) { const b = B[i]; const dx = b.x - px0, dy = b.y - py0, dz = b.z - pz0; d2[i] = dx * dx + dy * dy + dz * dz; }
    ord.sort(function (i, j) { return d2[j] - d2[i]; });
    const maxWH = Math.max(W, H);
    let grads = 0, mode = 'lighter';
    for (let k = 0; k < ord.length && grads < 300; k++) {
      const b = B[ord[k]];
      const dx = b.x - px0, dy = b.y - py0, dz = b.z - pz0;
      const vz = dx * fx + dy * fy + dz * fz;
      if (vz <= 0) continue;
      const sr = (b.r / vz) * focal;
      if (sr < 0.7 || sr > 3 * maxWH) continue;
      const sx = cx + (dx * rx + dy * ry + dz * rz) * focal / vz;
      const sy = cy - (dx * ux + dy * uy + dz * uz) * focal / vz;
      if (sx + sr < 0 || sx - sr > W || sy + sr < 0 || sy - sr > H) continue;
      const near = clamp(Math.sqrt(d2[ord[k]]) / b.r - 0.6, 0, 1);     /* inside a cloud you do not see it */
      const al = b.a * detail * near;
      if (al < 0.002) continue;
      const want = b.dust ? 'source-over' : 'lighter';
      if (want !== mode) { mode = want; ctx.globalCompositeOperation = want; }
      if (sr < 2.5) { ctx.fillStyle = 'rgba(' + b.col + ',' + (al * 0.9).toFixed(3) + ')'; ctx.fillRect(sx - sr, sy - sr, sr * 2, sr * 2); continue; }
      grads++;
      const g = ctx.createRadialGradient(sx, sy, 0, sx, sy, sr);
      g.addColorStop(0, 'rgba(' + b.col + ',' + al.toFixed(3) + ')');
      g.addColorStop(0.4, 'rgba(' + b.col + ',' + (al * 0.5).toFixed(3) + ')');
      g.addColorStop(1, 'rgba(' + b.col + ',0)');
      ctx.fillStyle = g;
      ctx.fillRect(sx - sr, sy - sr, sr * 2, sr * 2);
    }
    /* stars */
    if (mode !== 'lighter') ctx.globalCompositeOperation = 'lighter';
    const pa = detail * clamp(galPx / 350, 0.15, 1);
    const size = galPx > 2500 ? 2 : galPx > 900 ? 1.5 : 1;
    const near2 = 1500 * LY_AU * 1500 * LY_AU;
    const N = MW.n, PX = MW.px, PY = MW.py, PZ = MW.pz, BS = MW.bucketStart;
    for (let b = 0; b < 8; b++) {
      const st = MW.styles[b];
      ctx.fillStyle = 'rgba(' + st.col + ',' + (st.a * pa).toFixed(3) + ')';
      for (let i = BS[b]; i < BS[b + 1]; i++) {
        const dx = PX[i] - px0, dy = PY[i] - py0, dz = PZ[i] - pz0;
        const vz = dx * fx + dy * fy + dz * fz;
        if (vz <= 0) continue;
        if (dx * dx + dy * dy + dz * dz < near2) continue;
        const sx = cx + (dx * rx + dy * ry + dz * rz) * focal / vz;
        if (sx < 0 || sx > W) continue;
        const sy = cy - (dx * ux + dy * uy + dz * uz) * focal / vz;
        if (sy < 0 || sy > H) continue;
        ctx.fillRect(sx, sy, size, size);
      }
    }
  }
  ctx.restore();
  /* the arms, once they are big enough to read */
  if (galPx > 150 && galPx < 6000 && detail > 0.4 && cam.fov > FOV_DEFAULT * 0.98) {
    ctx.save();
    ctx.font = '600 9.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    if ('letterSpacing' in ctx) ctx.letterSpacing = '1.5px';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    const selArm = selectedId && byId[selectedId] && byId[selectedId].armIndex;
    const selBar = selectedId === 'bar';
    for (let i = 0; i < MW.armLabels.length; i++) {
      const L = MW.armLabels[i];
      const q = project(L.pos); if (!q || q.x < 30 || q.x > W - 30 || q.y < 30 || q.y > H - 30) continue;
      const name = (L.alt || STRUCTURE.arms.find(function (s) { return s.arm === L.arm; }).name).toUpperCase();
      const hot = L.arm == null ? selBar : (selArm != null && selArm === L.arm);
      const tw = ctx.measureText(name).width;
      const box = { x: q.x - tw / 2 - 3, y: q.y - 7, w: tw + 6, h: 14 };
      if (overlapsLabel(box) && !hot) continue;
      backdropLabelBoxes.push(box);
      ctx.fillStyle = hot ? 'rgba(255,240,200,0.98)' : 'rgba(190,205,255,' + (0.6 * fade).toFixed(3) + ')';
      ctx.fillText(name, q.x, q.y);
    }
    ctx.restore();
  }
  /* a name, once the whole thing fits on screen */
  if (c && galPx > 2 && galPx < 2600) {
    const small = galPx < 24;
    const ly = small ? c.y : c.y - galPx * 0.62 - 6, lx = small ? c.x + galPx + 6 : c.x;
    if (lx > 40 && lx < W - 40 && ly > 30 && ly < H - 30) {
      ctx.save();
      ctx.font = '600 11px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
      const tw = ctx.measureText('Milky Way').width;
      const box = small ? { x: lx, y: ly - 8, w: tw + 8, h: 16 } : { x: lx - tw / 2 - 4, y: ly - 8, w: tw + 8, h: 16 };
      if (!overlapsLabel(box)) {
        backdropLabelBoxes.push(box);
        ctx.textAlign = small ? 'left' : 'center'; ctx.textBaseline = 'middle';
        ctx.fillStyle = 'rgba(220,215,245,' + (0.85 * fade).toFixed(3) + ')';
        ctx.fillText('Milky Way', lx, ly);
      }
      ctx.restore();
    }
  }
}

function drawSky() {
  const fade = skyFade();
  /* The bulk catalogue is real geometry, not painted backdrop, so it must not fade
     out with the procedural sky when the camera leaves the solar neighbourhood. */
  if (bulkStars && opts.stars) {
    /* from outside the galaxy the catalogue stars are a bright knot at the Sun; let them recede */
    const camR = Math.hypot(camPos.x, camPos.y, camPos.z) / LY_AU;
    drawBulkStars(1 - 0.85 * clamp((camR - 15000) / 40000, 0, 1));
  }
  if (fade <= 0.01) return;

  if (opts.milkyway) {
    ctx.save();
    let mode = null;
    for (let i = 0; i < mwBlobs.length; i++) {
      const b = mwBlobs[i];
      const p = projectDir(b.d);
      if (!p) continue;
      const r = b.ang * focal;
      if (r > 4 * Math.max(W, H)) continue;            /* magnified sky: the glow is a flat wash */
      if (p.x < -r || p.x > W + r || p.y < -r || p.y > H + r) continue;
      const want = b.dust ? 'source-over' : 'lighter';
      if (want !== mode) { ctx.globalCompositeOperation = want; mode = want; }
      const a = b.a * fade * opts.mwGain * (deep.on ? deep.galaxyFade : 1);
      const soft = b.col.replace('rgb(', 'rgba(').replace(')', ',');
      const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, r);
      g.addColorStop(0, soft + a.toFixed(4) + ')');
      g.addColorStop(0.42, soft + (a * 0.46).toFixed(4) + ')');
      g.addColorStop(1, soft + '0)');
      ctx.fillStyle = g;
      ctx.fillRect(p.x - r, p.y - r, r * 2, r * 2);
    }
    ctx.restore();
  }
  if (!opts.stars) return;
  if (bulkStars) return;   /* already drawn above, unfaded */

  let bucket = -1;
  let streaks = null;
  for (let i = 0; i < skyStars.length; i++) {
    const s = skyStars[i];
    const p = projectDir(s.d);
    if (!p || p.x < 0 || p.x > W || p.y < 0 || p.y > H) { s.px = null; continue; }
    const b = Math.min(5, (s.m * 6) | 0);
    if (b !== bucket) {
      bucket = b;
      const a = (0.16 + b * 0.17) * fade;
      ctx.fillStyle = 'rgba(255,255,255,' + a.toFixed(3) + ')';
    }
    /* if it moved far since the last frame, trail it instead of dotting it */
    if (s.px != null) {
      const dx = p.x - s.px, dy = p.y - s.py;
      const d2 = dx * dx + dy * dy;
      if (d2 > 3 && d2 < 90000) {
        if (!streaks) { streaks = true; ctx.beginPath(); }
        ctx.moveTo(s.px, s.py);
        ctx.lineTo(p.x, p.y);
        s.px = p.x; s.py = p.y;
        continue;
      }
    }
    ctx.fillRect(p.x, p.y, s.s, s.s);
    s.px = p.x; s.py = p.y;
  }
  if (streaks) {
    ctx.strokeStyle = 'rgba(255,255,255,' + (0.34 * fade).toFixed(3) + ')';
    ctx.lineWidth = 1;
    ctx.stroke();
  }
}

/* 40,000 real stars. The array is sorted brightest-first at build time, so
   magnitude buckets change monotonically and fillStyle is set ~6 times a
   frame rather than 40,000 times. */
function drawBulkStars(fade) {
  const S = bulkStars, n = S.n;
  const px = camPos.x, py = camPos.y, pz = camPos.z;
  const fx = fwd.x, fy = fwd.y, fz = fwd.z;
  const rx = right.x, ry = right.y, rz = right.z;
  const ux = up.x, uy = up.y, uz = up.z;
  /* within the solar neighbourhood the catalogue magnitude is right as it
     stands; only once you actually travel does it need recomputing */
  const camR = Math.hypot(px, py, pz);
  const nearSun = camR < 3000;
  let bucket = -1, drawn = 0;
  const budget = 30000;
  const nDraw = opts.allStars ? n : S.defaultN;
  const life = deep.on ? starLives() : null, dY = deep.Y;

  for (let i = 0; i < nDraw && drawn < budget; i++) {
    if (life && life[i] < dY) continue;                 /* this star has already died */
    const dx = S.x[i] - px, dy = S.y[i] - py, dz = S.z[i] - pz;
    const vz = dx * fx + dy * fy + dz * fz;
    if (vz <= 1e-6) continue;
    const sx = cx + (dx * rx + dy * ry + dz * rz) * focal / vz;
    if (sx < 0 || sx > W) continue;
    const sy = cy - (dx * ux + dy * uy + dz * uz) * focal / vz;
    if (sy < 0 || sy > H) continue;

    let m = S.mag[i] * 0.01;
    if (!nearSun) {
      const dpc = Math.sqrt(dx * dx + dy * dy + dz * dz) / PC_AU;
      const dcat = Math.sqrt(S.x[i] * S.x[i] + S.y[i] * S.y[i] + S.z[i] * S.z[i]) / PC_AU;
      m += 5 * (Math.log(dpc / Math.max(dcat, 1e-6)) / Math.LN10);
      if (m > 9) continue;
    }
    drawn++;

    const b = m < 1.2 ? 0 : m < 2.6 ? 1 : m < 4.0 ? 2 : m < 5.2 ? 3 : m < 6.4 ? 4 : 5;
    if (b >= 2) {
      if (b !== bucket) {
        bucket = b;
        const a = [0, 0, 0.82, 0.52, 0.30, 0.16][b] * fade;
        ctx.fillStyle = 'rgba(255,255,255,' + a.toFixed(3) + ')';
      }
      ctx.fillRect(sx, sy, b <= 3 ? 1.5 : 1, b <= 3 ? 1.5 : 1);
    } else {
      /* the few hundred bright ones get their real colour and a little bloom */
      bucket = -1;
      const col = ciColor(S.ci[i] * 0.001);
      const r = b === 0 ? 2.3 : 1.6;
      const gr = b === 0 ? 11 : 6;
      const g = ctx.createRadialGradient(sx, sy, 0, sx, sy, gr);
      g.addColorStop(0, hexA(col, (0.55 * fade).toFixed(3)));
      g.addColorStop(1, hexA(col, 0));
      ctx.fillStyle = g;
      ctx.fillRect(sx - gr, sy - gr, gr * 2, gr * 2);
      ctx.fillStyle = col;
      ctx.beginPath(); ctx.arc(sx, sy, r, 0, TAU); ctx.fill();
    }
  }
}

/* 12,000 NGC/IC objects as faint smudges on the celestial sphere */
function drawBulkDso() {
  if (!bulkDso || !opts.ngc || (deep.on && deep.Y > 1e8)) return;
  const fade = skyFade();
  if (fade <= 0.03) return;
  const D = bulkDso, n = D.n, INV = 1 / 32767;
  const GRAD_BUDGET = 400;
  let gradients = 0;
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  for (let i = 0; i < n; i++) {
    const dxr = D.x[i] * INV, dyr = D.y[i] * INV, dzr = D.z[i] * INV;
    const vz = dxr * fwd.x + dyr * fwd.y + dzr * fwd.z;
    if (vz <= 0.02) continue;
    const sx = cx + (dxr * right.x + dyr * right.y + dzr * right.z) * focal / vz;
    if (sx < -20 || sx > W + 20) continue;
    const sy = cy - (dxr * up.x + dyr * up.y + dzr * up.z) * focal / vz;
    if (sy < -20 || sy > H + 20) continue;
    const mag = D.mag[i] * 0.01;
    const a = clamp((15.5 - mag) / 22, 0.05, 0.5) * fade * 0.55;
    /* MajAx is the major axis (a diameter), so halve it for a radius */
    const rr = clamp((D.maj[i] * 0.01 / 120) * DEG * focal, 1.2, 90 * FOV_DEFAULT / cam.fov);   /* true size once magnified */
    const t = D.typeNames[D.type[i]];
    const col = (t === 'G' || t === 'GPair' || t === 'GTrpl' || t === 'GGroup') ? '#cdd6ee'
              : (t === 'GCl') ? '#ffe6c4' : (t === 'OCl') ? '#dbe6ff'
              : (t === 'PN') ? '#8fd8c8' : '#d89a9a';
    /* A gradient per object is far too expensive at 12,000 rows. Only the ones
       big enough for the falloff to be visible get one; the rest are flat dots. */
    if (rr < 2.6 || gradients >= GRAD_BUDGET) {
      ctx.fillStyle = hexA(col, a.toFixed(3));
      ctx.fillRect(sx - rr, sy - rr, rr * 2, rr * 2);
    } else {
      gradients++;
      const g = ctx.createRadialGradient(sx, sy, 0, sx, sy, rr);
      g.addColorStop(0, hexA(col, a.toFixed(3)));
      g.addColorStop(0.5, hexA(col, (a * 0.35).toFixed(3)));
      g.addColorStop(1, hexA(col, 0));
      ctx.fillStyle = g;
      ctx.fillRect(sx - rr, sy - rr, rr * 2, rr * 2);
    }
  }
  ctx.restore();
}

function drawBelts() {
  if (deep.on && deep.Y > 1e15) return;
  const days = simJd - J2000;
  for (let f = 0; f < beltFields.length; f++) {
    const field = beltFields[f];
    const cfg = field.cfg;
    if (cfg.id === 'asteroidbelt' && bulkSbdb) continue;   /* real JPL asteroids replace the procedural belt */
    if (cfg.spherical && !opts.oort) continue;
    if (!cfg.spherical && !opts.belts) continue;
    if (cfg.maxCamAU && cam.dist > cfg.maxCamAU) continue;
    if (cfg.minCamAU && cam.dist < cfg.minCamAU) continue;

    const col = cfg.color;
    const rgb = [parseInt(col.slice(1, 3), 16), parseInt(col.slice(3, 5), 16), parseInt(col.slice(5, 7), 16)];
    let bucket = -1;
    for (let k = 0; k < field.n; k++) {
      const th = field.th[k] + field.w[k] * days;
      const r = field.r[k];
      let px, py, pz;
      if (cfg.spherical) {
        px = field.ux[k] * r; py = field.uy[k] * r; pz = field.uz[k] * r;
      } else {
        const c = Math.cos(th) * r, s = Math.sin(th) * r;
        px = c * field.ux[k] + s * field.vx[k];
        py = c * field.uy[k] + s * field.vy[k];
        pz = c * field.uz[k] + s * field.vz[k];
      }
      /* projection inlined: project() allocates an argument object and a result
         object, which at ~6,600 particles a frame is pure garbage-collector churn */
      const ddx = px - camPos.x, ddy = py - camPos.y, ddz = pz - camPos.z;
      const vz = ddx * fwd.x + ddy * fwd.y + ddz * fwd.z;
      if (vz <= 1e-11) continue;
      const sxp = cx + (ddx * right.x + ddy * right.y + ddz * right.z) * focal / vz;
      if (sxp < 0 || sxp > W) continue;
      const syp = cy - (ddx * up.x + ddy * up.y + ddz * up.z) * focal / vz;
      if (syp < 0 || syp > H) continue;
      const b = Math.min(3, (field.b[k] * 4) | 0);
      if (b !== bucket) {
        bucket = b;
        const a = 0.09 + b * 0.10;
        ctx.fillStyle = 'rgba(' + rgb[0] + ',' + rgb[1] + ',' + rgb[2] + ',' + a.toFixed(3) + ')';
      }
      ctx.fillRect(sxp, syp, 1, 1);
    }
  }
}

/* orbit ellipse, sampled directly in eccentric anomaly (no solving needed) */
function drawOrbit(o, T, alpha) {
  const el = o.el;
  const a = el.a + (el.aD || 0) * T, e = el.e + (el.eD || 0) * T;
  const I = (el.i + (el.iD || 0) * T) * DEG;
  const w = (el.w + (el.wD || 0) * T) * DEG;
  const O = (el.O + (el.OD || 0) * T) * DEG;
  const argP = w - O;
  const N = 220;
  ctx.beginPath();
  let started = false, drawn = 0;
  for (let k = 0; k <= N; k++) {
    const p3 = orbitalToEcl(a, e, (k / N) * TAU, argP, O, I);
    const p = project(p3);
    if (!p) { started = false; continue; }
    if (Math.abs(p.x) > 60000 || Math.abs(p.y) > 60000) { started = false; continue; }
    if (!started) { ctx.moveTo(p.x, p.y); started = true; } else { ctx.lineTo(p.x, p.y); }
    drawn++;
  }
  if (!drawn) return;
  ctx.strokeStyle = o.color;
  ctx.globalAlpha = alpha;
  ctx.lineWidth = 1;
  ctx.stroke();
  ctx.globalAlpha = 1;
}

function drawMoonOrbit(o, alpha) {
  const par = byId[o.parent];
  const p3 = o.refPlane === 'ecliptic' ? { x: 0, y: 0, z: 1 } : par.pole;
  let u = cross({ x: 0, y: 0, z: 1 }, p3);
  if (Math.hypot(u.x, u.y, u.z) < 1e-6) u = { x: 1, y: 0, z: 0 };
  u = norm(u);
  const v = cross(p3, u);
  const ci = Math.cos(o.inc), si = Math.sin(o.inc);
  const N = 96;
  ctx.beginPath();
  let started = false;
  for (let k = 0; k <= N; k++) {
    const th = (k / N) * TAU;
    const c = Math.cos(th) * o.aAU, s = Math.sin(th) * o.aAU;
    const p = project({
      x: par.pos.x + c * u.x + s * (v.x * ci + p3.x * si),
      y: par.pos.y + c * u.y + s * (v.y * ci + p3.y * si),
      z: par.pos.z + c * u.z + s * (v.z * ci + p3.z * si)
    });
    if (!p) { started = false; continue; }
    if (!started) { ctx.moveTo(p.x, p.y); started = true; } else { ctx.lineTo(p.x, p.y); }
  }
  ctx.strokeStyle = o.color;
  ctx.globalAlpha = alpha;
  ctx.lineWidth = 1;
  ctx.stroke();
  ctx.globalAlpha = 1;
}

/* a lit sphere: day side toward the Sun, dark limb away from it */
/* ---- Body frame: where a body's pole and prime meridian point right now ----
   IAU convention: the prime meridian is measured from the node Q of the body's
   equator on the ICRF equator (Q is 90° east of the pole's right ascension), and
   turns by W = W0 + 360/rotationDays · days since J2000. Tidally locked moons keep
   their prime meridian toward the planet instead. */
function bodyFrame(o) {
  const pole = o.pole || { x: 0, y: 0, z: 1 };
  let m;
  if (o.tidal && o.parent && byId[o.parent]) {
    const q = byId[o.parent].pos;
    m = { x: q.x - o.pos.x, y: q.y - o.pos.y, z: q.z - o.pos.z };
  } else {
    const d = o.data || {};
    const ra = (d.poleRA != null ? d.poleRA : 0) * DEG;
    const Q = eqToEcl(-Math.sin(ra), Math.cos(ra), 0);
    const W = ((d.W0 || 0) + (d.rotationDays ? 360 / d.rotationDays : 0) * (simJd - J2000)) * DEG;
    const pq = cross(pole, Q);
    m = { x: Q.x * Math.cos(W) + pq.x * Math.sin(W), y: Q.y * Math.cos(W) + pq.y * Math.sin(W), z: Q.z * Math.cos(W) + pq.z * Math.sin(W) };
  }
  const dd = m.x * pole.x + m.y * pole.y + m.z * pole.z;
  m = norm({ x: m.x - dd * pole.x, y: m.y - dd * pole.y, z: m.z - dd * pole.z });
  return { pole: pole, m: m, e: cross(pole, m) };
}
/* a point on the body's surface (lat/lon in degrees) -> unit vector in the world */
function surfVec(F, lat, lon) {
  const cl = Math.cos(lat * DEG), sl = Math.sin(lat * DEG), cn = Math.cos(lon * DEG), sn = Math.sin(lon * DEG);
  return { x: cl * (cn * F.m.x + sn * F.e.x) + sl * F.pole.x, y: cl * (cn * F.m.y + sn * F.e.y) + sl * F.pole.y, z: cl * (cn * F.m.z + sn * F.e.z) + sl * F.pole.z };
}
/* the drawing frame for a disc: two axes orthogonal to the line of sight, so the
   visible hemisphere's edge lands exactly on the drawn circle */
function discFrame(toCam) {
  const dd = right.x * toCam.x + right.y * toCam.y + right.z * toCam.z;
  let e1 = { x: right.x - dd * toCam.x, y: right.y - dd * toCam.y, z: right.z - dd * toCam.z };
  if (Math.hypot(e1.x, e1.y, e1.z) < 1e-6) { const du = up.x * toCam.x + up.y * toCam.y + up.z * toCam.z; e1 = { x: up.x - du * toCam.x, y: up.y - du * toCam.y, z: up.z - du * toCam.z }; }
  e1 = norm(e1);
  return { toCam: toCam, e1: e1, e2: cross(toCam, e1) };
}
function subLon(F, dir) { return Math.atan2(dir.x * F.e.x + dir.y * F.e.y + dir.z * F.e.z, dir.x * F.m.x + dir.y * F.m.y + dir.z * F.m.z) / DEG; }
/* Draw a closed polygon of unit vectors on the sphere, clipped to the visible
   hemisphere. Hidden stretches are replaced by the limb arc, so nothing folds
   across the disc. Polygons must wind counter-clockwise as seen from outside. */
function spherePoly(p, r, DF, verts) {
  const n = verts.length, tc = DF.toCam;
  let anyVis = false, allVis = true;
  const dep = new Array(n);
  for (let i = 0; i < n; i++) { const v = verts[i]; dep[i] = v.x * tc.x + v.y * tc.y + v.z * tc.z; if (dep[i] > 0) anyVis = true; else allVis = false; }
  if (!anyVis) return false;
  const px = function (v) { return { x: p.x + (v.x * DF.e1.x + v.y * DF.e1.y + v.z * DF.e1.z) * r, y: p.y - (v.x * DF.e2.x + v.y * DF.e2.y + v.z * DF.e2.z) * r }; };
  ctx.beginPath();
  if (allVis) {
    for (let i = 0; i < n; i++) { const q = px(verts[i]); if (i === 0) ctx.moveTo(q.x, q.y); else ctx.lineTo(q.x, q.y); }
    ctx.closePath(); return true;
  }
  let started = false, pendingExit = null, firstEntry = null;
  for (let i = 0; i < n; i++) {
    const A = verts[i], B = verts[(i + 1) % n], da = dep[i], db = dep[(i + 1) % n];
    if (da > 0) { const q = px(A); if (started) ctx.lineTo(q.x, q.y); else { ctx.moveTo(q.x, q.y); started = true; } }
    if ((da > 0) !== (db > 0)) {
      const t = da / (da - db);
      const c = norm({ x: A.x + t * (B.x - A.x), y: A.y + t * (B.y - A.y), z: A.z + t * (B.z - A.z) });
      const q = px(c), ang = Math.atan2(q.y - p.y, q.x - p.x);
      if (da > 0) { ctx.lineTo(q.x, q.y); pendingExit = ang; }
      else if (pendingExit != null) { ctx.arc(p.x, p.y, r, pendingExit, ang, true); pendingExit = null; }
      else { if (started) ctx.lineTo(q.x, q.y); else { ctx.moveTo(q.x, q.y); started = true; } firstEntry = ang; }
    }
  }
  if (pendingExit != null && firstEntry != null) ctx.arc(p.x, p.y, r, pendingExit, firstEntry, true);
  ctx.closePath();
  return true;
}
/* orient a polygon counter-clockwise as seen from outside the sphere */
function orientCCW(verts) {
  let cx = 0, cy = 0, cz = 0;
  verts.forEach(function (v) { cx += v.x; cy += v.y; cz += v.z; });
  const c = norm({ x: cx, y: cy, z: cz });
  let s = 0;
  for (let i = 0; i < verts.length; i++) { const A = verts[i], B = verts[(i + 1) % verts.length]; const x = cross(A, B); s += x.x * c.x + x.y * c.y + x.z * c.z; }
  return s < 0 ? verts.slice().reverse() : verts;
}
const bandCache = {};
/* the band between two latitudes as a counter-clockwise polygon (east along the south edge, west along the north edge) */
function bandVerts(F, lat0, lat1) {
  const N = 48, out = [];
  if (lat0 > -89.9) for (let k = 0; k <= N; k++) out.push(surfVec(F, lat0, -180 + 360 * k / N));
  else out.push(surfVec(F, -90, 0));
  if (lat1 < 89.9) for (let k = N; k >= 0; k--) out.push(surfVec(F, lat1, -180 + 360 * k / N));
  else out.push(surfVec(F, 90, 0));
  return out;
}
function spotVerts(F, sp, lon) {
  const N = 28, out = [], stretch = 1 / Math.max(0.15, Math.cos(sp.lat * DEG));
  for (let k = 0; k < N; k++) { const t = k / N * TAU; out.push(surfVec(F, sp.lat + Math.sin(t) * sp.h / 2, (lon != null ? lon : sp.lon) + Math.cos(t) * sp.w / 2 * stretch)); }
  return out;
}
function polyVerts(F, pts) { return orientCCW(pts.map(function (q) { return surfVec(F, q[0], q[1]); })); }
const craterRng = mulberry32(77);
const CRATERS = []; for (let k = 0; k < 40; k++) CRATERS.push({ lat: (craterRng() * 2 - 1) * 75, lon: craterRng() * 360, s: 1.5 + craterRng() * 5 });
const cloudRng = mulberry32(31);
const CLOUDS = []; for (let k = 0; k < 40; k++) CLOUDS.push({ lat: (cloudRng() * 2 - 1) * 65, lon: cloudRng() * 360, w: 14 + cloudRng() * 30, h: 4 + cloudRng() * 7, a: 0.35 + cloudRng() * 0.4 });

/* ---- A lit body: base colour, its surface, then day and night ---- */
function paintSphere(p, r, o, sunScreen) {
  const host = o.kind === 'exoplanet' ? (byId[o.lightId] || byId[o.parent]) : (o.lightId && byId[o.lightId]);
  const lightPos = host ? host.pos : { x: 0, y: 0, z: 0 };
  const L = norm({ x: lightPos.x - o.pos.x, y: lightPos.y - o.pos.y, z: lightPos.z - o.pos.z });
  const toCam = norm({ x: camPos.x - o.pos.x, y: camPos.y - o.pos.y, z: camPos.z - o.pos.z });
  const DF = discFrame(toCam);
  const Lx = L.x * DF.e1.x + L.y * DF.e1.y + L.z * DF.e1.z;
  const Ly = -(L.x * DF.e2.x + L.y * DF.e2.y + L.z * DF.e2.z);
  const Lz = L.x * toCam.x + L.y * toCam.y + L.z * toCam.z;      /* > 0: the Sun is behind us */
  const S = o.dry ? null : (SURFACES[o.id] || null);
  const F = (S && opts.surfaces && r > 8) ? bodyFrame(o) : null;

  ctx.save();
  ctx.beginPath(); ctx.arc(p.x, p.y, r, 0, TAU); ctx.clip();
  ctx.fillStyle = (S && S.ocean) || (S && S.base) || o.color;
  ctx.fillRect(p.x - r, p.y - r, r * 2, r * 2);
  if (!F) {                                                   /* small, or the simple look: a shaded ball */
    const g = ctx.createRadialGradient(p.x + Lx * r * 0.45, p.y + Ly * r * 0.45, r * 0.04, p.x, p.y, r * 1.02);
    g.addColorStop(0, o.color); g.addColorStop(0.55, o.color2 || o.color); g.addColorStop(1, 'rgba(0,0,0,0.92)');
    ctx.fillStyle = g; ctx.fillRect(p.x - r, p.y - r, r * 2, r * 2);
    if (o.bands && r > 7) {                                   /* soft cloud bands for the giants */
      const polePt = project({ x: o.pos.x + o.pole.x * o.radiusAU, y: o.pos.y + o.pole.y * o.radiusAU, z: o.pos.z + o.pole.z * o.radiusAU });
      const ang = polePt ? Math.atan2(polePt.y - p.y, polePt.x - p.x) + Math.PI / 2 : 0;
      ctx.save(); ctx.translate(p.x, p.y); ctx.rotate(ang);
      for (let b = -4; b <= 4; b++) { const y = (b / 5) * r, h = r * 0.13; ctx.fillStyle = (b % 2 === 0) ? 'rgba(255,255,255,0.07)' : 'rgba(0,0,0,0.12)'; ctx.fillRect(-r, y - h / 2, r * 2, h); }
      ctx.restore();
    }
  } else {
    if (S.bands) S.bands.forEach(function (b) { if (spherePoly(p, r, DF, bandVerts(F, b[0], b[1]))) { ctx.fillStyle = b[2]; ctx.fill(); } });
    if (S.land) {
      ctx.fillStyle = S.landColor;
      S.land.forEach(function (poly) { if (spherePoly(p, r, DF, polyVerts(F, poly))) ctx.fill(); });
      ctx.fillStyle = S.desert;
      S.deserts.forEach(function (sp) { if (spherePoly(p, r, DF, spotVerts(F, sp))) ctx.fill(); });
      ctx.fillStyle = S.ice;
      if (spherePoly(p, r, DF, polyVerts(F, S.greenlandIce))) ctx.fill();
    }
    if (S.iceCaps) { ctx.fillStyle = S.ice; S.iceCaps.forEach(function (c) { if (spherePoly(p, r, DF, bandVerts(F, c[0], c[1]))) ctx.fill(); }); }
    if (S.spots) S.spots.forEach(function (sp) {
      let lon = sp.lon;
      if (sp.sys2Lon != null) {                             /* a feature that keeps still in System II: convert to the System III frame we rotate */
        const west3 = sp.sys2Lon + (sp.driftDegDay || 0) * (simJd - sp.epochJd) + 0.266 * (simJd - J2000) + 241.65;
        lon = -west3;
      }
      if (spherePoly(p, r, DF, spotVerts(F, sp, lon))) { ctx.fillStyle = sp.color; ctx.fill(); }
    });
    if (S.craters && r > 20) {
      for (let k = 0; k < CRATERS.length; k++) {
        const c = CRATERS[k], v = surfVec(F, c.lat, c.lon);
        if (v.x * toCam.x + v.y * toCam.y + v.z * toCam.z <= 0.05) continue;
        const qx = p.x + (v.x * DF.e1.x + v.y * DF.e1.y + v.z * DF.e1.z) * r, qy = p.y - (v.x * DF.e2.x + v.y * DF.e2.y + v.z * DF.e2.z) * r;
        const cs = c.s * r / 100;
        ctx.fillStyle = 'rgba(0,0,0,0.22)'; ctx.beginPath(); ctx.arc(qx, qy, cs, 0, TAU); ctx.fill();
        ctx.fillStyle = 'rgba(255,255,255,0.10)'; ctx.beginPath(); ctx.arc(qx - cs * 0.3, qy - cs * 0.3, cs * 0.6, 0, TAU); ctx.fill();
      }
    }
    if (S.clouds && r > 14) {
      ctx.fillStyle = 'rgba(255,255,255,0.5)';
      for (let k = 0; k < S.clouds; k++) { const c = CLOUDS[k]; if (spherePoly(p, r, DF, spotVerts(F, c))) { ctx.globalAlpha = c.a * 0.6; ctx.fill(); } }
      ctx.globalAlpha = 1;
    }
    if (o.rings && r > 24) ringShadowOnPlanet(o, p, r, F, DF, L);
  }
  /* day and night: the terminator is a great circle, an ellipse on screen */
  const m = Math.hypot(Lx, Ly);
  ctx.beginPath();
  ctx.arc(p.x, p.y, r, 0, TAU);
  if (m < 1e-4) {
    if (Lz > 0) ctx.rect(p.x - r, p.y - r, r * 2, r * 2);     /* Sun straight behind us: fully lit, cancel the disc */
  } else {
    litPath(p, r, Math.atan2(Ly, Lx), Lz);
  }
  ctx.fillStyle = 'rgba(0,0,0,0.9)';
  ctx.fill('evenodd');
  /* shading across the lit side: brightest under the Sun, dim at the terminator */
  const sx = p.x + Lx * r * 0.9, sy = p.y + Ly * r * 0.9;
  const g = ctx.createRadialGradient(sx, sy, 0, sx, sy, r * 1.9);
  g.addColorStop(0, 'rgba(255,255,255,0.10)'); g.addColorStop(0.35, 'rgba(0,0,0,0)'); g.addColorStop(1, 'rgba(0,0,0,0.78)');
  ctx.fillStyle = g; ctx.fillRect(p.x - r, p.y - r, r * 2, r * 2);
  ctx.restore();
  /* a thin atmosphere on the lit limb */
  const atm = S && opts.surfaces && S.atmosphere;
  if (atm && r > 12) {
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    const ga = ctx.createRadialGradient(p.x, p.y, r * 0.96, p.x, p.y, r * 1.08);
    ga.addColorStop(0, hexA(atm, 0.0)); ga.addColorStop(0.35, hexA(atm, 0.28)); ga.addColorStop(1, hexA(atm, 0));
    ctx.fillStyle = ga;
    if (m >= 1e-4) { ctx.beginPath(); litPath(p, r * 1.08, Math.atan2(Ly, Lx), Lz); ctx.clip(); ctx.fillRect(p.x - r * 1.1, p.y - r * 1.1, r * 2.2, r * 2.2); }
    else if (Lz > 0) { ctx.beginPath(); ctx.arc(p.x, p.y, r * 1.08, 0, TAU); ctx.fill(); }
    ctx.restore();
  }
}
/* the sunlit region of a disc: the half toward the Sun plus (or minus) the terminator ellipse */
function litPath(p, r, ang, Lz) {
  ctx.moveTo(p.x + Math.cos(ang - Math.PI / 2) * r, p.y + Math.sin(ang - Math.PI / 2) * r);
  ctx.arc(p.x, p.y, r, ang - Math.PI / 2, ang + Math.PI / 2);
  ctx.ellipse(p.x, p.y, Math.max(0.01, r * Math.abs(Lz)), r, ang, Math.PI / 2, 3 * Math.PI / 2, Lz < 0);
  ctx.closePath();
}
function ringZones(rg) {
  return rg.gapInnerKm
    ? [[rg.innerKm, 92000, 0.16, 0.20], [92000, rg.gapInnerKm, 0.62, 0.62], [rg.gapOuterKm, rg.outerKm, 0.42, 0.45]]
    : [[rg.innerKm, rg.outerKm, 0.30, 0.30]];
}
function ringShadowOnPlanet(o, p, r, F, DF, L) {
  const rg = o.rings, R = o.radiusAU, pole = F.pole;
  const lp = L.x * pole.x + L.y * pole.y + L.z * pole.z;
  if (Math.abs(lp) < 1e-4) return;
  const lonS = subLon(F, L);
  const zones = ringZones(rg);
  /* where each surface point's ray toward the Sun crosses the ring plane */
  const hitRadius = function (lat, lon) {
    const n = surfVec(F, lat, lon);
    const np = n.x * pole.x + n.y * pole.y + n.z * pole.z;
    const t = -np * R / lp;
    if (t <= 0) return -1;
    return Math.hypot(n.x * R + L.x * t, n.y * R + L.y * t, n.z * R + L.z * t);
  };
  for (let z = 0; z < zones.length; z++) {
    const inner = zones[z][0] / AU_KM, outer = zones[z][1] / AU_KM;
    const runs = [];                                     /* per longitude: [latA, latB] of the shadowed run */
    for (let k = 0; k <= 36; k++) {
      const lon = lonS - 90 + 180 * k / 36;
      let a = null, b = null;
      for (let lat = -85; lat <= 85; lat += 1) {
        const rho = hitRadius(lat, lon);
        const inZone = rho > inner && rho < outer;
        if (inZone) { if (a == null) a = lat; b = lat + 1; }
      }
      runs.push(a == null ? null : [lon, a, b]);
    }
    /* group consecutive longitudes into polygons */
    let seg = [];
    const flush = function () {
      if (seg.length >= 2) {
        const verts = [];
        seg.forEach(function (s) { verts.push(surfVec(F, s[1], s[0])); });
        for (let i = seg.length - 1; i >= 0; i--) verts.push(surfVec(F, seg[i][2], seg[i][0]));
        if (spherePoly(p, r, DF, verts)) { ctx.fillStyle = 'rgba(0,0,0,' + zones[z][3].toFixed(2) + ')'; ctx.fill(); }
      }
      seg = [];
    };
    runs.forEach(function (s) { if (s) seg.push(s); else flush(); });
    flush();
  }
}

/* ---- Rings: real annuli with the Cassini division, the planet's shadow on them,
   and the near half drawn over the disc. ---- */
function drawRings(o, p, planetR, sunScreen, front) {
  const rg = o.rings;
  const zones = ringZones(rg);
  const pole = o.pole;
  /* a far planet is drawn as a minimum-size dot; scale the rings with it so they never hide inside */
  const truePx = (o.radiusAU / p.z) * focal, scale = planetR > truePx ? planetR / truePx : 1;
  let u = cross({ x: 0, y: 0, z: 1 }, pole);
  if (Math.hypot(u.x, u.y, u.z) < 1e-6) u = { x: 1, y: 0, z: 0 };
  u = norm(u);
  const v = cross(pole, u);
  const lightPos = (o.lightId && byId[o.lightId]) ? byId[o.lightId].pos : { x: 0, y: 0, z: 0 };
  const L = norm({ x: lightPos.x - o.pos.x, y: lightPos.y - o.pos.y, z: lightPos.z - o.pos.z });
  const R = o.radiusAU, N = 96;
  const lp = L.x * pole.x + L.y * pole.y + L.z * pole.z;
  ctx.save();
  for (let z = 0; z < zones.length; z++) {
    const inner = zones[z][0] / AU_KM * scale, outer = zones[z][1] / AU_KM * scale, alpha = zones[z][2];
    const pts = [];
    for (let k = 0; k <= N; k++) {
      const th = (k / N) * TAU;
      const c = Math.cos(th), s = Math.sin(th);
      const dx = c * u.x + s * v.x, dy = c * u.y + s * v.y, dz = c * u.z + s * v.z;
      const pi = project({ x: o.pos.x + dx * inner, y: o.pos.y + dy * inner, z: o.pos.z + dz * inner });
      const po = project({ x: o.pos.x + dx * outer, y: o.pos.y + dy * outer, z: o.pos.z + dz * outer });
      if (!pi || !po) { pts.push(null); continue; }
      /* in the planet's shadow: the ring point's ray toward the Sun passes through the planet */
      const rm = (inner + outer) / 2, qx = dx * rm, qy = dy * rm, qz = dz * rm;
      const along = qx * L.x + qy * L.y + qz * L.z;
      const perp2 = qx * qx + qy * qy + qz * qz - along * along;
      const shadow = along < 0 && perp2 < R * R;
      pts.push({ i: pi, o: po, f: pi.z < p.z, shadow: shadow });
    }
    for (let k = 0; k < N; k++) {
      const A = pts[k], B = pts[k + 1];
      if (!A || !B || A.f !== front) continue;
      ctx.beginPath();
      ctx.moveTo(A.i.x, A.i.y); ctx.lineTo(A.o.x, A.o.y);
      ctx.lineTo(B.o.x, B.o.y); ctx.lineTo(B.i.x, B.i.y);
      ctx.closePath();
      /* the rings are lit from the Sun's side of the plane; from the other side they are dimmer */
      const cp = (camPos.x - o.pos.x) * pole.x + (camPos.y - o.pos.y) * pole.y + (camPos.z - o.pos.z) * pole.z;
      const sameSide = (cp > 0) === (lp > 0);
      ctx.fillStyle = hexA(rg.color, A.shadow || B.shadow ? alpha * 0.12 : alpha * (sameSide ? 1 : 0.55));
      ctx.fill();
    }
  }
  ctx.restore();
}

function drawBlackHole(o, p) {
  const shadowPx = Math.max(1.2, (o.radiusAU / p.z) * focal);
  const rot = (simJd - J2000) * 0.9;

  if (shadowPx < 3) {
    const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, 16);
    g.addColorStop(0, 'rgba(255,190,120,0.85)');
    g.addColorStop(0.35, 'rgba(255,140,60,0.28)');
    g.addColorStop(1, 'rgba(255,120,40,0)');
    ctx.fillStyle = g;
    ctx.fillRect(p.x - 16, p.y - 16, 32, 32);
    ctx.fillStyle = '#000';
    ctx.beginPath(); ctx.arc(p.x, p.y, Math.max(1, shadowPx), 0, TAU); ctx.fill();
    return;
  }
  /* accretion disc: a real 3D annulus so it tilts as the camera orbits */
  const pole = { x: 0.3, y: 0.2, z: 0.93 };
  let u = norm(cross({ x: 0, y: 0, z: 1 }, pole));
  if (!isFinite(u.x)) u = { x: 1, y: 0, z: 0 };
  const v = cross(pole, u);
  const inner = o.radiusAU * 2.3, outer = o.radiusAU * 9.5;
  const N = 120;

  /* sample the annulus once, keeping depth so we can split it around
     the shadow: far half, then the black hole, then the near half */
  const seg = [];
  for (let k = 0; k <= N; k++) {
    const th = (k / N) * TAU + rot;
    const c = Math.cos(th), s = Math.sin(th);
    const d = { x: c * u.x + s * v.x, y: c * u.y + s * v.y, z: c * u.z + s * v.z };
    const pi = project({ x: o.pos.x + d.x * inner, y: o.pos.y + d.y * inner, z: o.pos.z + d.z * inner });
    const po = project({ x: o.pos.x + d.x * outer, y: o.pos.y + d.y * outer, z: o.pos.z + d.z * outer });
    /* Doppler beaming: material sweeping toward the camera is brighter */
    const vel = { x: -s * u.x + c * v.x, y: -s * u.y + c * v.y, z: -s * u.z + c * v.z };
    const toward = -(vel.x * fwd.x + vel.y * fwd.y + vel.z * fwd.z);
    seg.push(pi && po ? { i: pi, o: po, th: th, beam: clamp(0.18 + 0.82 * (0.5 + 0.5 * toward), 0.12, 1), far: pi.z > p.z } : null);
  }

  function paintHalf(far) {
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    for (let k = 0; k < N; k++) {
      const A = seg[k], B = seg[k + 1];
      if (!A || !B || A.far !== far) continue;
      /* swirling filaments, so the disc reads as turbulent gas */
      const fil = 0.72 + 0.28 * Math.sin(A.th * 7.3 + rot * 2.1) * Math.sin(A.th * 3.1 - rot * 1.3);
      const b = ((A.beam + B.beam) / 2) * fil;
      const mx0 = (A.i.x + B.i.x) / 2, my0 = (A.i.y + B.i.y) / 2;
      const mx1 = (A.o.x + B.o.x) / 2, my1 = (A.o.y + B.o.y) / 2;
      const g = ctx.createLinearGradient(mx0, my0, mx1, my1);
      /* white-hot at the inner edge, cooling and thinning outward */
      g.addColorStop(0, 'rgba(255,246,226,' + (0.75 * b).toFixed(3) + ')');
      g.addColorStop(0.22, 'rgba(255,200,120,' + (0.50 * b).toFixed(3) + ')');
      g.addColorStop(0.60, 'rgba(255,140,55,' + (0.20 * b).toFixed(3) + ')');
      g.addColorStop(1, 'rgba(200,70,20,0)');
      ctx.beginPath();
      ctx.moveTo(A.i.x, A.i.y); ctx.lineTo(A.o.x, A.o.y);
      ctx.lineTo(B.o.x, B.o.y); ctx.lineTo(B.i.x, B.i.y);
      ctx.closePath();
      ctx.fillStyle = g;
      ctx.fill();
    }
    ctx.restore();
  }

  paintHalf(true);                       /* the far side of the disc */

  /* event horizon: a genuinely black hole punched through it */
  ctx.globalCompositeOperation = 'source-over';
  ctx.fillStyle = '#000';
  ctx.beginPath(); ctx.arc(p.x, p.y, shadowPx, 0, TAU); ctx.fill();

  /* photon ring hugging the shadow */
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  ctx.strokeStyle = 'rgba(255,214,164,0.9)';
  ctx.lineWidth = Math.max(1, shadowPx * 0.06);
  ctx.beginPath(); ctx.arc(p.x, p.y, shadowPx * 1.04, 0, TAU); ctx.stroke();
  ctx.restore();

  /* Gravitational lensing: light from the far side of the disc bends
     around the hole and shows up as arcs above and below the shadow,
     strongest perpendicular to the disc's projected major axis. */
  const a0 = seg[0], aMid = seg[(N / 2) | 0];
  let majAng = 0;
  if (a0 && aMid) majAng = Math.atan2(aMid.i.y - a0.i.y, aMid.i.x - a0.i.x);
  const r1 = shadowPx * 1.02, r2 = shadowPx * 2.25, M = 64;
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  const lg = ctx.createRadialGradient(p.x, p.y, r1, p.x, p.y, r2);
  lg.addColorStop(0, 'rgba(255,238,210,0.62)');
  lg.addColorStop(0.30, 'rgba(255,193,116,0.30)');
  lg.addColorStop(1, 'rgba(255,140,60,0)');
  ctx.fillStyle = lg;
  for (let k = 0; k < M; k++) {
    const t0 = (k / M) * TAU, t1 = ((k + 1) / M) * TAU;
    const wgt = Math.pow(Math.abs(Math.sin(t0 - majAng)), 1.7);
    if (wgt < 0.03) continue;
    ctx.globalAlpha = wgt;
    ctx.beginPath();
    ctx.moveTo(p.x + Math.cos(t0) * r1, p.y + Math.sin(t0) * r1);
    ctx.lineTo(p.x + Math.cos(t0) * r2, p.y + Math.sin(t0) * r2);
    ctx.lineTo(p.x + Math.cos(t1) * r2, p.y + Math.sin(t1) * r2);
    ctx.lineTo(p.x + Math.cos(t1) * r1, p.y + Math.sin(t1) * r1);
    ctx.closePath();
    ctx.fill();
  }
  ctx.restore();

  paintHalf(false);                      /* the near side, in front of the shadow */

  glow(p, shadowPx * 6, 'rgb(255,170,90)', 0.20);
}

/* A neutron star is ~22 km across, never resolvable in practice, so it is drawn as
   an intense point. Pulsars get the lighthouse beams they are known for, swept at a
   viewable rate rather than the real one — 716 Hz would just be a blur. */
function drawPulsarBeams(o, p) {
  /* the full beams belong to the pulsar you are looking at; from anywhere else a
     pulsar shows only a short flicker so the sky is not slashed by cones */
  const focused = selectedId === o.id || cam.follow === o.id || focusCategory === o.category;
  const maxL = focused ? 40000 : 34;
  const th = (simJd - J2000) * 1.7;
  const tilt = 0.55;                  /* magnetic axis offset from the spin axis */
  const ct = Math.cos(tilt), stl = Math.sin(tilt);
  const c = Math.cos(th), s = Math.sin(th);
  const m = norm({
    x: o.spinAxis.x * ct + (o.spinU.x * c + o.spinV.x * s) * stl,
    y: o.spinAxis.y * ct + (o.spinU.y * c + o.spinV.y * s) * stl,
    z: o.spinAxis.z * ct + (o.spinU.z * c + o.spinV.z * s) * stl
  });
  const len = Math.max(o.radiusAU * 80, p.z * 0.40);
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  for (let side = -1; side <= 1; side += 2) {
    const tp = project({ x: o.pos.x + m.x * len * side,
                         y: o.pos.y + m.y * len * side,
                         z: o.pos.z + m.z * len * side });
    if (!tp) continue;
    let dx = tp.x - p.x, dy = tp.y - p.y;
    let L = Math.hypot(dx, dy);
    if (!isFinite(L) || L < 2 || L > 40000) continue;
    if (L > maxL) { const k = maxL / L; dx *= k; dy *= k; tp.x = p.x + dx; tp.y = p.y + dy; L = maxL; }
    const nx = -dy / L, ny = dx / L;
    const w0 = 1.5, w1 = Math.min(L * 0.17, focused ? 240 : 5);
    const col = o.sub === 'magnetar' ? '210,160,255' : '150,205,255';
    const g = ctx.createLinearGradient(p.x, p.y, tp.x, tp.y);
    g.addColorStop(0, 'rgba(' + col + ',0.42)');
    g.addColorStop(0.35, 'rgba(' + col + ',0.13)');
    g.addColorStop(1, 'rgba(' + col + ',0)');
    ctx.fillStyle = g;
    ctx.beginPath();
    ctx.moveTo(p.x + nx * w0, p.y + ny * w0);
    ctx.lineTo(tp.x + nx * w1, tp.y + ny * w1);
    ctx.lineTo(tp.x - nx * w1, tp.y - ny * w1);
    ctx.lineTo(p.x - nx * w0, p.y - ny * w0);
    ctx.closePath();
    ctx.fill();
  }
  ctx.restore();
}

function drawRemnant(o, p) {
  const rPx = (o.radiusAU / p.z) * focal;
  const core = Math.max(rPx, opts.realSize ? 0.6 : o.minPx);
  const magnetar = o.sub === 'magnetar';
  if (o.beams) drawPulsarBeams(o, p);
  glow(p, Math.max(core * 7, 20), magnetar ? 'rgb(198,140,255)' : 'rgb(160,206,255)', 0.22);
  glow(p, Math.max(core * 2.6, 12), magnetar ? 'rgb(236,212,255)' : 'rgb(224,240,255)', 0.60);
  const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, core);
  g.addColorStop(0, '#ffffff');
  g.addColorStop(0.6, magnetar ? '#e6ccff' : '#e2f0ff');
  g.addColorStop(1, o.color);
  ctx.fillStyle = g;
  ctx.beginPath(); ctx.arc(p.x, p.y, core, 0, TAU); ctx.fill();
}

function drawDeepSky(o, p) {
  const rPx = Math.min(4200, Math.max(1.6, (o.radiusAU / p.z) * focal));
  const sub = o.sub;
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';

  if (sub === 'globular' || sub === 'open') {
    const halo = rPx * (sub === 'globular' ? 1.5 : 1.2);
    const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, halo);
    g.addColorStop(0, hexA(o.color, sub === 'globular' ? 0.34 : 0.24));
    g.addColorStop(0.45, hexA(o.color, 0.09));
    g.addColorStop(1, hexA(o.color, 0));
    ctx.fillStyle = g;
    ctx.fillRect(p.x - halo, p.y - halo, halo * 2, halo * 2);
    if (rPx > 2.2 && o.members) {
      const sz = clamp(rPx * 0.022, 1.1, 3.4);
      for (let k = 0; k < o.members.length; k++) {
        const m = o.members[k];
        const mx = p.x + m.x * rPx, my = p.y + m.y * rPx;
        if (mx < -20 || mx > W + 20 || my < -20 || my > H + 20) continue;
        if (m.b > 0.82 && sz > 1.6) {
          const gr = sz * 4;
          const gg = ctx.createRadialGradient(mx, my, 0, mx, my, gr);
          gg.addColorStop(0, hexA(o.color, 0.5));
          gg.addColorStop(1, hexA(o.color, 0));
          ctx.fillStyle = gg;
          ctx.fillRect(mx - gr, my - gr, gr * 2, gr * 2);
        }
        ctx.globalAlpha = m.b;
        ctx.fillStyle = '#ffffff';
        ctx.beginPath(); ctx.arc(mx, my, sz * (0.5 + m.b * 0.5), 0, TAU); ctx.fill();
        ctx.globalAlpha = 1;
      }
    }

  } else if (sub === 'planetary') {
    /* a shell seen from outside: bright rim, thinner middle */
    const g = ctx.createRadialGradient(p.x, p.y, rPx * 0.15, p.x, p.y, rPx);
    g.addColorStop(0, hexA(o.color2, 0.06));
    g.addColorStop(0.42, hexA(o.color2, 0.20));
    g.addColorStop(0.68, hexA(o.color, 0.42));
    g.addColorStop(0.92, hexA(o.color, 0.16));
    g.addColorStop(1, hexA(o.color, 0));
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(p.x, p.y, rPx, 0, TAU); ctx.fill();
    ctx.fillStyle = 'rgba(255,255,255,0.9)';
    ctx.beginPath(); ctx.arc(p.x, p.y, Math.max(0.8, rPx * 0.05), 0, TAU); ctx.fill();

  } else if (sub === 'galaxy' || sub === 'spiral') {
    const flat = o.data.flat != null ? o.data.flat : 0.4;
    ctx.translate(p.x, p.y);
    ctx.rotate(o.data.tilt || 0);
    ctx.scale(1, flat);
    const g = ctx.createRadialGradient(0, 0, 0, 0, 0, rPx);
    g.addColorStop(0, 'rgba(255,255,255,0.50)');
    g.addColorStop(0.14, hexA(o.color, 0.34));
    g.addColorStop(0.45, hexA(o.color, 0.13));
    g.addColorStop(1, hexA(o.color, 0));
    ctx.fillStyle = g;
    ctx.beginPath(); ctx.arc(0, 0, rPx, 0, TAU); ctx.fill();

  } else {
    for (let k = 0; k < o.puffs.length; k++) {
      const q = o.puffs[k];
      const rr = rPx * q.r;
      const px = p.x + q.x * rPx, py = p.y + q.y * rPx;
      const g = ctx.createRadialGradient(px, py, 0, px, py, rr);
      g.addColorStop(0, hexA(q.c, q.a));
      g.addColorStop(0.32, hexA(q.c, q.a * 0.62));
      g.addColorStop(0.62, hexA(q.c, q.a * 0.22));
      g.addColorStop(1, hexA(q.c, 0));
      ctx.fillStyle = g;
      ctx.fillRect(px - rr, py - rr, rr * 2, rr * 2);
    }
  }
  ctx.restore();
}

function glow(p, radius, color, alpha) {
  const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, radius);
  const soft = color.replace('rgb(', 'rgba(').replace(')', ',');
  g.addColorStop(0, color);
  g.addColorStop(0.16, soft + (alpha * 0.55).toFixed(3) + ')');
  g.addColorStop(0.38, soft + (alpha * 0.22).toFixed(3) + ')');
  g.addColorStop(0.65, soft + (alpha * 0.07).toFixed(3) + ')');
  g.addColorStop(1, 'rgba(255,200,120,0)');
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  ctx.globalAlpha = alpha;
  ctx.fillStyle = g;
  ctx.fillRect(p.x - radius, p.y - radius, radius * 2, radius * 2);
  ctx.restore();
}

function render(dt) {
  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, W, H);

  backdropLabelBoxes = []; labelBoxes = [];
  drawSky();
  drawBulkGalaxies();
  drawBulkDso();
  drawMilkyWay();
  drawRegions();
  drawBelts();
  drawSmallBodies(project({ x: 0, y: 0, z: 0 }));

  const T = (simJd - J2000) / 36525;
  const sunScreen = project({ x: 0, y: 0, z: 0 });

  /* orbit paths */
  if (opts.orbits) {
    for (let i = 0; i < objects.length; i++) {
      const o = objects[i];
      if (!o.el || o.kind === 'probe') continue;
      if (o.kind === 'dwarf' && cam.dist > 30000) continue;
      if (o.kind === 'dwarf' && o.category !== 'dwarfs' && selectedId !== o.id && cam.follow !== o.id) continue;
      let a = o.kind === 'planet' ? 0.30 : 0.18;
      if (deep.on) { if (deepGone(o) || deep.shot === 'earth') continue; if (o.deepBlur) a = 0.62; }   /* a planet moving too fast to see becomes a ring of light */
      drawOrbit(o, T, selectedId === o.id ? 0.75 : a);
    }
    for (let i = 0; i < objects.length; i++) {
      const o = objects[i];
      if (o.kind !== 'moon' || !opts.moons) continue;
      const par = byId[o.parent];
      const dp = Math.hypot(par.pos.x - camPos.x, par.pos.y - camPos.y, par.pos.z - camPos.z);
      if ((o.aAU / dp) * focal < 14) continue;
      drawMoonOrbit(o, 0.22);
    }
    for (let i = 0; i < objects.length; i++) {
      const o = objects[i];
      if (o.kind !== 'exoplanet' || !o.parent) continue;
      const host = byId[o.parent];
      if (!host) continue;
      const dp = Math.hypot(host.pos.x - camPos.x, host.pos.y - camPos.y, host.pos.z - camPos.z);
      if ((o.aAU / dp) * focal < 14) continue;
      const N = 96;
      ctx.beginPath();
      let started = false;
      for (let k = 0; k <= N; k++) {
        const th = (k / N) * TAU;
        const c = Math.cos(th) * o.aAU, sn = Math.sin(th) * o.aAU;
        const q = project({ x: host.pos.x + c * o.ou.x + sn * o.ov.x,
                            y: host.pos.y + c * o.ou.y + sn * o.ov.y,
                            z: host.pos.z + c * o.ou.z + sn * o.ov.z });
        if (!q) { started = false; continue; }
        if (!started) { ctx.moveTo(q.x, q.y); started = true; } else { ctx.lineTo(q.x, q.y); }
      }
      ctx.strokeStyle = o.color;
      ctx.globalAlpha = selectedId === o.id ? 0.6 : 0.24;
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    /* Parker Solar Probe's own ellipse */
    const parker = byId.parker;
    if (parker && parker.el && opts.probes && cam.dist < 60) drawOrbit(parker, T, 0.28);
  }

  /* L2 observatories: the path out from Earth, and the halo once parked */
  if (opts.probes && cam.dist < 40 && telescopeFocus()) {
    ctx.save();
    for (let i = 0; i < objects.length; i++) {
      const o = objects[i];
      if (o.mode !== 'l2halo' || o.launchJd == null) continue;
      const N = 90;
      ctx.beginPath();
      let started = false;
      for (let k = 0; k <= N; k++) {
        const t = o.launchJd + (o.arriveJd - o.launchJd) * (k / N);
        const q = l2Position(o, t);
        const sp = project(q);
        if (!sp) { started = false; continue; }
        if (!started) { ctx.moveTo(sp.x, sp.y); started = true; } else { ctx.lineTo(sp.x, sp.y); }
      }
      ctx.strokeStyle = o.color;
      ctx.globalAlpha = selectedId === o.id ? 0.85 : 0.35;
      ctx.setLineDash([4, 4]);
      ctx.lineWidth = 1.2;
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
    }
    ctx.restore();
  }

  drawRoute();

  /* spacecraft cruise lines */
  if (opts.probes && sunScreen) {
    ctx.save();
    for (let i = 0; i < objects.length; i++) {
      const o = objects[i];
      if (o.kind !== 'probe' || o.mode !== 'escape') continue;
      const p = project(o.pos);
      if (!p) continue;
      ctx.beginPath();
      ctx.moveTo(sunScreen.x, sunScreen.y);
      ctx.lineTo(p.x, p.y);
      ctx.strokeStyle = o.color;
      ctx.globalAlpha = selectedId === o.id ? 0.42 : 0.14;
      ctx.setLineDash([3, 5]);
      ctx.lineWidth = 1;
      ctx.stroke();
    }
    ctx.restore();
  }

  /* collect drawables. Clear every stale .screen first: a value left over from an
     earlier frame would put the selection reticle and pick() at a phantom position. */
  for (let i = 0; i < objects.length; i++) objects[i].screen = null;
  const camFromSunLy = Math.hypot(camPos.x, camPos.y, camPos.z) / LY_AU;
  const draw = [];
  for (let i = 0; i < objects.length; i++) {
    const o = objects[i];
    if (o.noRender) continue;
    if (o.kind === 'moon' && !opts.moons) continue;
    if (o.kind === 'probe' && !opts.probes) continue;
    if (o.kind === 'distant' && !deepVisible(o, camFromSunLy)) continue;
    if (deep.on && deepHidden(o)) continue;
    /* faint dwarf galaxies would only clutter the sky from home */
    if (o.data && o.data.dwarf && camFromSunLy < 40000 && selectedId !== o.id && focusCategory !== 'galaxies' && selectedId !== 'localgroup') continue;
    const p = project(o.pos);
    if (!p) continue;
    if (p.x < -200 || p.x > W + 200 || p.y < -200 || p.y > H + 200) continue;
    o.screen = p;
    draw.push(o);
  }
  draw.sort((a, b) => b.screen.z - a.screen.z);

  labelBoxes = backdropLabelBoxes.slice();
  const labels = [];

  for (let i = 0; i < draw.length; i++) {
    const o = draw[i];
    const p = o.screen;

    const farOut = camFromSunLy > 20000 && selectedId !== o.id && cam.follow !== o.id;
    if (o.kind === 'distant') { drawDistant(o, p); labels.push(o); continue; }
    if (o.kind === 'blackhole') {
      if (farOut && o.id !== 'sgrA') continue;
      drawBlackHole(o, p);
      labels.push(o);
      continue;
    }

    if (o.kind === 'deepsky') {
      drawDeepSky(o, p);
      labels.push(o);
      continue;
    }

    if (o.kind === 'neutron') {
      if (farOut) continue;
      drawRemnant(o, p);
      labels.push(o);
      continue;
    }

    if (o.kind === 'star' || o.kind === 'sun') {
      /* brightness as seen from where the camera actually is */
      let rPx = (o.radiusAU / p.z) * focal;
      /* brightness follows true range, not the view-axis depth p.z — otherwise a
         star appears to brighten simply for drifting toward the frame edge */
      const rng = Math.hypot(o.pos.x - camPos.x, o.pos.y - camPos.y, o.pos.z - camPos.z);
      const dpc = rng / (LY_AU * PC_LY);
      const mag = (o.kind === 'sun' ? 4.83 : o.absMag) + 5 * Math.log10(Math.max(dpc, 1e-12)) - 5;
      /* how much brighter than the naked-eye limit, compressed hard so a
         nearby star glows without washing out everything beside it */
      if (o.kind === 'star' && mag > 11 && farOut) continue;      /* too faint to see from here */
      const excess = Math.max(0, 6.5 - mag);
      const bright = clamp(excess / 9, 0.05, 1);
      const core = Math.max(opts.realSize ? 0.7 : (o.kind === 'sun' && camFromSunLy > 60000 ? 1.6 : o.minPx), rPx);
      const glowR = Math.max(core * 2.2, clamp(3.5 + 5.2 * Math.sqrt(excess), 4, 40));

      if (rPx > 2.0) {
        /* a resolved star: wide atmosphere, corona, tight halo, then the disc */
        const gk = (deep.on && o.kind === 'sun') ? 0.12 : 1;     /* a giant Sun would wash out the screen */
        glow(p, Math.max(rPx * 34 * gk, 170), 'rgb(255,214,160)', clamp(0.06 + bright * 0.14, 0.06, 0.20));
        glow(p, Math.max(rPx * 13 * gk, 70), 'rgb(255,230,185)', clamp(0.20 + bright * 0.40, 0.20, 0.60));
        glow(p, Math.max(rPx * 4.6, 26), 'rgb(255,248,228)', clamp(0.35 + bright * 0.50, 0.35, 0.85));
        const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, rPx);
        g.addColorStop(0, '#ffffff');
        g.addColorStop(0.55, (deep.on && o.kind === 'sun') ? o.color : '#fffaf0');
        g.addColorStop(0.88, o.color);
        g.addColorStop(1, o.kind === 'sun' ? '#ffab3d' : o.color);
        ctx.fillStyle = g;
        ctx.beginPath(); ctx.arc(p.x, p.y, rPx, 0, TAU); ctx.fill();
      } else {
        glow(p, glowR, 'rgb(255,230,190)', clamp(0.10 + bright * 0.32, 0.07, 0.42));
        ctx.fillStyle = '#fff';
        ctx.beginPath(); ctx.arc(p.x, p.y, core, 0, TAU); ctx.fill();
        ctx.fillStyle = o.color;
        ctx.globalAlpha = 0.75;
        ctx.beginPath(); ctx.arc(p.x, p.y, core * 1.7, 0, TAU); ctx.fill();
        ctx.globalAlpha = 1;
      }
      if (o.kind === 'sun' || mag < 5.2 || selectedId === o.id) labels.push(o);
      continue;
    }

    /* solid bodies */
    const litFrom = (o.lightId && byId[o.lightId]) ? project(byId[o.lightId].pos) : sunScreen;
    let rPx = (o.radiusAU / p.z) * focal;
    const r = opts.realSize ? Math.max(rPx, 0.6) : Math.max(rPx, o.minPx);

    if (o.kind === 'probe') {
      ctx.fillStyle = o.color;
      ctx.beginPath(); ctx.arc(p.x, p.y, 2.6, 0, TAU); ctx.fill();
      ctx.strokeStyle = o.color;
      ctx.globalAlpha = 0.55;
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(p.x, p.y, 6.5, 0, TAU); ctx.stroke();
      ctx.globalAlpha = 1;
      labels.push(o);
      continue;
    }

    if (o.rings && rPx > 1.2) drawRings(o, p, r, litFrom, false);

    if (r > 1.9) {
      paintSphere(p, r, o, litFrom);
    } else {
      ctx.fillStyle = o.color;
      ctx.beginPath(); ctx.arc(p.x, p.y, r, 0, TAU); ctx.fill();
      ctx.globalAlpha = 0.35;
      ctx.beginPath(); ctx.arc(p.x, p.y, r * 2.6, 0, TAU); ctx.fill();
      ctx.globalAlpha = 1;
    }

    if (o.rings && rPx > 1.2) drawRings(o, p, r, litFrom, true);

    if (o.kind === 'planet' || o.kind === 'dwarf' || o.kind === 'exoplanet') labels.push(o);
    else if (o.kind === 'moon' && r > 1.6) labels.push(o);
  }

  /* selection reticle */
  if (selectedId && byId[selectedId] && byId[selectedId].screen) {
    const o = byId[selectedId], p = o.screen;
    const rr = Math.max(11, Math.min(60, (o.radiusAU / p.z) * focal * 1.9));
    ctx.strokeStyle = 'rgba(160,175,255,0.9)';
    ctx.lineWidth = 1.2;
    for (let q = 0; q < 4; q++) {
      const a0 = q * Math.PI / 2 + Math.PI / 8;
      ctx.beginPath();
      ctx.arc(p.x, p.y, rr, a0, a0 + Math.PI / 4);
      ctx.stroke();
    }
  }

  drawSatellites();
  if (opts.labels) drawLabels(labels);
  drawDeepFields();
  drawConstellationNames();
}

/* Which labels earn a place at the current zoom. Everything is still there and
   searchable — but a label should mean "this matters here", not "this exists". */
function labelWanted(o) {
  if (selectedId === o.id || cam.follow === o.id) return true;
  const p = o.screen; if (!p) return false;
  /* from outside the galaxy everything near the Sun shares one pixel; the galaxy gets the name */
  if (o.kind !== 'deepsky' && o.kind !== 'ngc' && o.kind !== 'constellation' && o.kind !== 'deepfield' && o.kind !== 'distant' &&
      Math.hypot(camPos.x, camPos.y, camPos.z) > 4e5 * LY_AU) return false;
  const rPx = (o.radiusAU / p.z) * focal;
  switch (o.kind) {
    case 'sun': return Math.hypot(camPos.x, camPos.y, camPos.z) < 4e5 * LY_AU;   /* from outside, the galaxy gets the name */
    case 'planet': return true;
    case 'dwarf': return o.category === 'dwarfs' || rPx > 2.5;
    case 'moon': return rPx > 1.6;
    case 'exoplanet': return rPx > 1.5;
    case 'star': return (o.appMag != null ? o.appMag : 9) < (cam.dist > 1000 ? 2.5 : 1.2) || rPx > 2;
    case 'deepsky': case 'ngc': return focusCategory === o.category || (rPx > 6 && cam.dist > 3000) || (o.lg && selectedId === 'localgroup');
    case 'probe':
      if (o.parent) return selectedId === o.parent || cam.follow === o.parent;
      if (o.mode === 'l2' || o.mode === 'l2halo' || o.mode === 'trailing') return cam.dist < 0.6 || focusCategory === 'telescopes';
      return cam.dist < 300;
    case 'distant': return focusCategory === 'deepuniverse' || cam.dist > 2e7 * LY_AU || o.sub === 'cluster' && cam.dist > 2e6 * LY_AU;
    case 'neutron': case 'blackhole':
      if (Math.hypot(camPos.x, camPos.y, camPos.z) > 4e5 * LY_AU) return false;   /* from outside, the galaxy gets the name */
      return focusCategory === o.category || cam.dist < 3 || o.id === 'sgrA';
    default: return true;
  }
}
function drawLabels(list) {
  list = list.filter(labelWanted);
  list.sort((a, b) => (b.prio || 0) - (a.prio || 0));
  const maxLabels = cam.dist > 200 ? 16 : cam.dist > 5 ? 22 : 30;
  ctx.font = '11.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  ctx.textBaseline = 'middle';
  let count = 0;
  for (let i = 0; i < list.length && count < maxLabels; i++) {
    const o = list[i], p = o.screen;
    if (!p || p.x < 4 || p.x > W - 4 || p.y < 4 || p.y > H - 4) continue;
    const rPx = Math.max(o.minPx || 2, (o.radiusAU / p.z) * focal);
    const off = Math.min(26, rPx + 7);
    const tw = ctx.measureText(o.name).width;
    const box = { x: p.x + off, y: p.y - 7, w: tw + 6, h: 14 };
    if (overlapsLabel(box) && selectedId !== o.id) continue;
    labelBoxes.push(box);
    count++;
    const sel = selectedId === o.id;
    ctx.fillStyle = sel ? 'rgba(190,200,255,1)' : 'rgba(226,226,236,0.82)';
    ctx.fillText(o.name, box.x, p.y);
    ctx.strokeStyle = sel ? 'rgba(160,175,255,0.6)' : 'rgba(190,190,210,0.25)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(p.x + off - 5, p.y); ctx.lineTo(p.x + Math.min(off, rPx + 2), p.y);
    ctx.stroke();
  }
}

/* Deep fields and telescope trajectories stay out of the sky unless you are
   actually looking at telescopes — selecting one, or browsing the category. */
function telescopeFocus() {
  const sel = selectedId && byId[selectedId];
  if (cam.fov < FOV_DEFAULT * 0.98) return { scope: null };
  if (sel && sel.kind === 'deepfield') return { scope: null };
  if (sel && sel.category === 'telescopes') return { scope: (sel.data && sel.data.scope) || null };
  if (focusCategory === 'telescopes' || focusCategory === 'deepfields') return { scope: null };
  return null;
}

function drawDeepFields() {
  if (!opts.deepfields || deep.on) return;
  const focus = telescopeFocus();
  if (!focus) return;
  const fade = skyFade();
  if (fade <= 0.03) return;
  const selObj = selectedId && byId[selectedId];
  ctx.save();
  ctx.font = '600 10px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  ctx.textBaseline = 'middle';
  ctx.textAlign = 'left';
  for (let i = 0; i < deepFieldList.length; i++) {
    const d = deepFieldList[i], f = d.data;
    /* when one telescope is selected, show only the patches it looked at */
    if (focus.scope && f.scopes.indexOf(focus.scope) < 0) continue;
    const p = projectDir(d.dir);
    if (!p) continue;
    /* the patch at its true angular size — a fraction of a pixel at normal zoom */
    const hw = Math.tan(f.arcmin / 120 * DEG) * focal;
    const hh = Math.tan((f.arcminH || f.arcmin) / 120 * DEG) * focal;
    const ext = Math.max(hw, hh);
    const big = ext > 12;
    if (!big && (p.x < 40 || p.x > W - 40 || p.y < 24 || p.y > H - 24)) continue;
    if (big && (p.x + ext < -20 || p.x - ext > W + 20 || p.y + ext < -20 || p.y - ext > H + 20)) continue;
    d.screen = { x: p.x, y: p.y, z: 1 };
    const sel = selObj === d;
    const col = sel ? 'rgba(200,220,255,0.95)' : hexA(f.color, (0.55 * fade).toFixed(3));
    let lx = p.x + 15, ly = p.y;
    if (!big) {
      /* corner brackets: the real field is far smaller than this marker */
      const r = 10, arm = 4;
      ctx.strokeStyle = col; ctx.lineWidth = 1.2;
      const corners = [[-1, -1], [1, -1], [1, 1], [-1, 1]];
      for (let c = 0; c < 4; c++) {
        const sx = corners[c][0], sy = corners[c][1];
        ctx.beginPath();
        ctx.moveTo(p.x + sx * r, p.y + sy * r - sy * arm);
        ctx.lineTo(p.x + sx * r, p.y + sy * r);
        ctx.lineTo(p.x + sx * r - sx * arm, p.y + sy * r);
        ctx.stroke();
      }
      if (sel) { drawMoonForScale(p, p.x - r, p.x + r); drawFieldCaption(f, p.x - r, p.y + r); }
    } else {
      const fr = skyFrame(f.ra, f.dec, p, f.pa || 0);
      const c = patchCorners(p, fr, hw, hh);
      if (sel && ext > 24) drawFieldContents(f, p, fr, hw, hh);
      ctx.strokeStyle = col; ctx.lineWidth = sel ? 1.4 : 1;
      ctx.setLineDash(f.planned ? [5, 4] : []);
      ctx.beginPath(); ctx.moveTo(c[0].x, c[0].y);
      for (let k = 1; k < 4; k++) ctx.lineTo(c[k].x, c[k].y);
      ctx.closePath(); ctx.stroke(); ctx.setLineDash([]);
      let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
      c.forEach(function (q) { minX = Math.min(minX, q.x); minY = Math.min(minY, q.y); maxX = Math.max(maxX, q.x); maxY = Math.max(maxY, q.y); });
      lx = minX; ly = minY - 10;
      if (sel) {
        drawFieldHighlights(f);
        drawMoonForScale(p, minX, maxX);
        drawFieldCaption(f, minX, maxY);
      }
    }
    ctx.font = '600 10px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const tw = ctx.measureText(d.name).width;
    const box = { x: lx, y: ly - 8, w: tw + 8, h: 16 };
    if (!overlapsLabel(box) || sel) {
      labelBoxes.push(box);
      ctx.fillStyle = sel ? 'rgba(210,226,255,0.98)' : hexA(f.color, (0.78 * fade).toFixed(3));
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText(d.name, lx, ly);
    }
  }
  ctx.restore();
}

/* Screen-space unit vectors for celestial east and north at a point on the sky,
   optionally rotated by a position angle (east of north) for tilted survey strips. */
function skyFrame(raH, decDeg, p, paDeg) {
  const dN = norm(raDecToVec(raH, decDeg + 0.05, 1));
  const dE = norm(raDecToVec(raH + 0.05 / 15 / Math.max(Math.cos(decDeg * DEG), 0.05), decDeg, 1));
  const pN = projectDir(dN), pE = projectDir(dE);
  let n = pN ? norm2(pN.x - p.x, pN.y - p.y) : { x: 0, y: -1 };
  let e = pE ? norm2(pE.x - p.x, pE.y - p.y) : { x: -1, y: 0 };
  if (paDeg) {
    const cs = Math.cos(paDeg * DEG), sn = Math.sin(paDeg * DEG);
    const n2 = { x: n.x * cs + e.x * sn, y: n.y * cs + e.y * sn };
    const e2 = { x: e.x * cs - n.x * sn, y: e.y * cs - n.y * sn };
    n = n2; e = e2;
  }
  return { e: e, n: n };
}
function norm2(x, y) { const m = Math.hypot(x, y) || 1; return { x: x / m, y: y / m }; }
function patchCorners(p, fr, hw, hh) {
  const out = [];
  [[-1, -1], [1, -1], [1, 1], [-1, 1]].forEach(function (s) {
    out.push({ x: p.x + fr.e.x * s[0] * hw + fr.n.x * s[1] * hh, y: p.y + fr.e.y * s[0] * hw + fr.n.y * s[1] * hh });
  });
  return out;
}

/* ---- What the telescope saw: an impression of the field at its real density ----
   We do not ship the photographs, so the patch is filled with synthetic galaxies —
   the right number for the field, sizes and colours drawn from what such images
   contain, clustered a little as real galaxies are. It is labelled as an
   impression on screen. Real highlights (GN-z11, Earendel…) sit at their true
   coordinates on top. */
function fieldSeed(str) { let h = 2166136261; for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 16777619); } return h >>> 0; }
function makeRng(seed) {
  let s = seed >>> 0;
  return function () { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}
function gaussian(r) { return Math.sqrt(-2 * Math.log(1 - r() + 1e-12)) * Math.cos(TAU * r()); }
const FIELD_COLS = [[168, 190, 255], [236, 240, 255], [246, 222, 170], [232, 140, 105]];
const FIELD_STYLES = FIELD_COLS.map(function (c) { return [0.28, 0.45, 0.65, 0.9].map(function (al) { return 'rgba(' + c.join(',') + ',' + al + ')'; }); });
function fieldContents(f) {
  if (f._gal) return f._gal;
  const r = makeRng(fieldSeed(f.id));
  const n = Math.min(f.count || 0, 20000);
  const webb = f.scopes.indexOf('JWST') >= 0;
  const clumps = [];
  for (let c = 0; c < 6; c++) clumps.push({ u: r() - 0.5, v: r() - 0.5, s: 0.03 + r() * 0.06 });
  const gal = new Float32Array(n * 7);
  for (let i = 0; i < n; i++) {
    let u, v;
    if (r() < 0.35) {
      const c = clumps[(r() * clumps.length) | 0];
      u = c.u + gaussian(r) * c.s; v = c.v + gaussian(r) * c.s;
      if (Math.abs(u) > 0.5 || Math.abs(v) > 0.5) { u = r() - 0.5; v = r() - 0.5; }
    } else { u = r() - 0.5; v = r() - 0.5; }
    const size = Math.exp(Math.log(0.5) + gaussian(r) * 0.75);          /* arcsec, log-normal */
    const t = r();
    const cls = t < (webb ? 0.30 : 0.42) ? 0 : t < 0.62 ? 1 : t < (webb ? 0.78 : 0.86) ? 2 : 3;
    gal[i * 7] = u; gal[i * 7 + 1] = v; gal[i * 7 + 2] = size;
    gal[i * 7 + 3] = 0.35 + r() * 0.65; gal[i * 7 + 4] = r() * Math.PI;
    gal[i * 7 + 5] = cls; gal[i * 7 + 6] = Math.pow(r(), 2.2);
  }
  const area = f.arcmin * (f.arcminH || f.arcmin);
  const nStars = n ? Math.min(40, Math.max(3, Math.round(area * 1.2))) : 0;
  const stars = [];
  for (let i = 0; i < nStars; i++) stars.push({ u: r() - 0.5, v: r() - 0.5, b: 0.3 + r() * 0.7 });
  f._gal = { gal: gal, n: n, stars: stars };
  return f._gal;
}
function drawFieldContents(f, p, fr, hw, hh) {
  const C = fieldContents(f);
  if (!C.n) return;
  const pxPerArcsec = (DEG / 3600) * focal;
  const c = patchCorners(p, fr, hw, hh);
  ctx.save();
  ctx.beginPath(); ctx.moveTo(c[0].x, c[0].y);
  for (let k = 1; k < 4; k++) ctx.lineTo(c[k].x, c[k].y);
  ctx.closePath(); ctx.clip();
  const g = C.gal, ex = fr.e.x * 2 * hw, ey = fr.e.y * 2 * hw, nx = fr.n.x * 2 * hh, ny = fr.n.y * 2 * hh;
  const base = Math.atan2(fr.e.y, fr.e.x);
  let cur = null, grads = 0;
  for (let i = 0; i < C.n; i++) {
    const o = i * 7;
    const x = p.x + ex * g[o] + nx * g[o + 1], y = p.y + ey * g[o] + ny * g[o + 1];
    if (x < -8 || x > W + 8 || y < -8 || y > H + 8) continue;
    const s = g[o + 2] * pxPerArcsec, cls = g[o + 5] | 0, b = g[o + 6];
    const lvl = b < 0.25 ? 0 : b < 0.5 ? 1 : b < 0.8 ? 2 : 3;
    const st = FIELD_STYLES[cls][lvl];
    if (st !== cur) { cur = st; ctx.fillStyle = st; }
    if (s < 1.3) { ctx.fillRect(x, y, 1, 1); continue; }
    if (s > 7 && grads < 60) {
      grads++;
      const col = FIELD_COLS[cls];
      const rg = ctx.createRadialGradient(x, y, 0, x, y, s / 2);
      rg.addColorStop(0, 'rgba(' + col.join(',') + ',' + (0.55 + b * 0.45).toFixed(2) + ')');
      rg.addColorStop(0.45, 'rgba(' + col.join(',') + ',' + (0.25 + b * 0.3).toFixed(2) + ')');
      rg.addColorStop(1, 'rgba(' + col.join(',') + ',0)');
      ctx.fillStyle = rg; cur = null;
    }
    ctx.beginPath();
    ctx.ellipse(x, y, s / 2, (s / 2) * g[o + 3], base + g[o + 4], 0, TAU);
    ctx.fill();
  }
  /* a few foreground Milky Way stars, with the telescope's own diffraction pattern */
  const spikes = f.scopes.indexOf('JWST') >= 0 ? 6 : f.scopes.indexOf('Hubble') >= 0 ? 4 : 0;
  for (let i = 0; i < C.stars.length; i++) {
    const s = C.stars[i];
    const x = p.x + ex * s.u + nx * s.v, y = p.y + ey * s.u + ny * s.v;
    if (x < -50 || x > W + 50 || y < -50 || y > H + 50) continue;
    const core = clamp((0.4 + s.b * 0.8) * pxPerArcsec, 1.5, 6);
    const len = clamp((4 + s.b * 14) * pxPerArcsec, 8, 170);
    if (spikes) {
      ctx.lineWidth = 1;
      const a0 = spikes === 6 ? Math.PI / 6 : Math.PI / 4;
      for (let k = 0; k < spikes; k++) {
        const ang = a0 + k * TAU / spikes, tx = x + Math.cos(ang) * len, ty = y + Math.sin(ang) * len;
        const lg = ctx.createLinearGradient(x, y, tx, ty);
        lg.addColorStop(0, 'rgba(255,245,225,' + (0.5 * s.b).toFixed(2) + ')'); lg.addColorStop(1, 'rgba(255,245,225,0)');
        ctx.strokeStyle = lg; ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(tx, ty); ctx.stroke();
      }
      if (spikes === 6) {   /* Webb's two short horizontal spikes from the secondary-mirror struts */
        for (let k = -1; k <= 1; k += 2) {
          const lg = ctx.createLinearGradient(x, y, x + k * len * 0.4, y);
          lg.addColorStop(0, 'rgba(255,245,225,' + (0.35 * s.b).toFixed(2) + ')'); lg.addColorStop(1, 'rgba(255,245,225,0)');
          ctx.strokeStyle = lg; ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + k * len * 0.4, y); ctx.stroke();
        }
      }
    }
    glow({ x: x, y: y }, core * 4, 'rgb(255,240,220)', 0.35 * s.b);
    ctx.fillStyle = '#fff8ee'; ctx.beginPath(); ctx.arc(x, y, core, 0, TAU); ctx.fill();
  }
  ctx.restore();
}
function drawFieldHighlights(f) {
  const hs = f.highlights || [];
  ctx.save();
  ctx.font = '600 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
  for (let i = 0; i < hs.length; i++) {
    const h = hs[i]; if (!h.dir) continue;
    const q = projectDir(h.dir);
    if (!q || q.x < -20 || q.x > W + 20 || q.y < -20 || q.y > H + 20) continue;
    const hot = f._shown === i;
    ctx.strokeStyle = hot ? '#ffe9a8' : 'rgba(255,210,122,0.9)'; ctx.lineWidth = hot ? 1.6 : 1.1;
    ctx.beginPath(); ctx.arc(q.x, q.y, hot ? 9 : 6, 0, TAU); ctx.stroke();
    if (hot) { ctx.beginPath(); ctx.arc(q.x, q.y, 14, 0, TAU); ctx.strokeStyle = 'rgba(255,233,168,0.35)'; ctx.stroke(); }
    const txt = h.name + (h.z != null ? ' · z ' + h.z : '');
    const tw = ctx.measureText(txt).width;
    const box = { x: q.x + 11, y: q.y - 8, w: tw + 8, h: 16 };
    if (!overlapsLabel(box) || hot) {
      labelBoxes.push(box);
      ctx.fillStyle = 'rgba(255,225,160,0.95)';
      ctx.fillText(txt, q.x + 11, q.y);
    }
  }
  ctx.restore();
}
function drawMoonForScale(p, minX, maxX) {
  const r = Math.tan(15.5 / 60 * DEG) * focal;      /* the full Moon is 31 arcmin across */
  if (r < 18 || r > 0.42 * Math.min(W, H)) return;
  let cxm = maxX + 26 + r;
  if (cxm + r > W - 10) cxm = minX - 26 - r;
  if (cxm - r < 10) return;
  ctx.save();
  ctx.strokeStyle = 'rgba(255,255,255,0.35)'; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.arc(cxm, p.y, r, 0, TAU); ctx.stroke(); ctx.setLineDash([]);
  ctx.font = '500 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  ctx.fillStyle = 'rgba(255,255,255,0.6)'; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  ctx.fillText('Full Moon, for scale', cxm, p.y + r + 6);
  ctx.restore();
}
function drawFieldCaption(f, minX, maxY) {
  ctx.save();
  ctx.textAlign = 'left'; ctx.textBaseline = 'top';
  ctx.font = '600 11px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  ctx.fillStyle = 'rgba(220,230,255,0.9)';
  const size = f.arcminH ? f.arcmin + '′ × ' + f.arcminH + '′' : (f.arcmin >= 90 ? (f.arcmin / 60).toFixed(1) + '° across' : f.arcmin + '′ across');
  ctx.fillText(f.name + ' · ' + size, minX, maxY + 10);
  if (f.count) {
    const drawn = Math.min(f.count, 20000);
    ctx.font = '500 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    ctx.fillStyle = 'rgba(200,210,235,0.62)';
    ctx.fillText((f.count >= 1000 ? f.count.toLocaleString() : 'thousands of') + ' galaxies' + (drawn < f.count ? ' (' + drawn.toLocaleString() + ' drawn)' : '') +
                 ' — an impression at the real density, not the photograph', minX, maxY + 25);
  } else if (f.planned) {
    ctx.font = '500 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    ctx.fillStyle = 'rgba(200,210,235,0.62)';
    ctx.fillText('Planned — nothing observed here yet', minX, maxY + 25);
  }
  ctx.restore();
}

/* turn to look along a direction, and set the magnification */
function lookAtDir(d, fovRad) {
  let dy = Math.atan2(-d.y, -d.x) - cam.yawGoal;
  dy -= TAU * Math.round(dy / TAU);
  cam.flight = null;
  cam.follow = null;      /* turning to look is not following */
  updateFocusChip();
  cam.yawGoal += dy;
  cam.pitchGoal = clamp(-Math.asin(clamp(d.z, -1, 1)), -1.45, 1.45);
  cam.rollGoal = 0; cam.pendYaw = cam.pendPitch = 0;          /* a telescope looks with a level horizon */
  cam.fovGoal = clamp(fovRad || FOV_DEFAULT, FOV_MIN, FOV_DEFAULT);
  idle = 0;
}
function fmtAngle(rad) {
  const am = rad / DEG * 60;
  if (am < 1) return Math.round(am * 60) + '″';
  if (am < 90) return am.toFixed(am < 10 ? 1 : 0) + '′';
  return (am / 60).toFixed(1) + '°';
}
function angSepArcmin(a, b) { return Math.acos(clamp(a.x * b.x + a.y * b.y + a.z * b.z, -1, 1)) / DEG * 60; }

/* ---- Cosmology: what a redshift means in years and light years ----
   Flat ΛCDM with the Planck 2018 numbers. Integrals are done in ln(1+z) so the
   very-high-z end stays well behaved. */
const COSMO = { H0: 67.4, Om: 0.315, Or: 9.2e-5 };
COSMO.OL = 1 - COSMO.Om - COSMO.Or;
COSMO.tH = 977.792 / COSMO.H0;                       /* Hubble time, billion years */
COSMO.dH = 299792.458 / COSMO.H0 * 3.261564e-3;      /* Hubble distance, billion light years */
function cosmoE(z) { const a = 1 + z; return Math.sqrt(COSMO.Or * a * a * a * a + COSMO.Om * a * a * a + COSMO.OL); }
function simpson(fn, a, b, n) {
  n = n & ~1; const h = (b - a) / n; let s = fn(a) + fn(b);
  for (let i = 1; i < n; i++) s += fn(a + i * h) * (i & 1 ? 4 : 2);
  return s * h / 3;
}
COSMO.age0 = COSMO.tH * simpson(function (x) { return 1 / cosmoE(Math.exp(x) - 1); }, 0, Math.log(1e6), 4000);
COSMO.horizon = COSMO.dH * simpson(function (x) { return Math.exp(x) / cosmoE(Math.exp(x) - 1); }, 0, Math.log(1e6), 4000);
function cosmo(z) {
  const X = Math.log(1 + z);
  const lb = simpson(function (x) { return 1 / cosmoE(Math.exp(x) - 1); }, 0, X, 400);
  const dc = simpson(function (x) { return Math.exp(x) / cosmoE(Math.exp(x) - 1); }, 0, X, 400);
  return { lookback: COSMO.tH * lb, age: COSMO.age0 - COSMO.tH * lb, comoving: COSMO.dH * dc };
}
function fmtGyr(g) {
  if (g >= 1) return g.toFixed(g < 10 ? 2 : 1) + ' billion years';
  if (g >= 0.001) return Math.round(g * 1000) + ' million years';
  return Math.round(g * 1e6).toLocaleString() + ' thousand years';
}
function fmtGly(g) { return g >= 1 ? g.toFixed(g < 10 ? 2 : 1) + ' billion light years' : Math.round(g * 1000) + ' million light years'; }
function cosmoSentence(z) {
  const c = cosmo(z);
  return 'Light left it ' + fmtGyr(c.lookback) + ' ago, when the universe was ' + fmtGyr(c.age) +
         ' old. Space has expanded since, so it is now about ' + fmtGly(c.comoving) + ' away.';
}


/* --- The deep universe: quasars, clusters, explosions, the first galaxies, and the
   walls and voids of the cosmic web. Built here, after the cosmology calculator,
   because redshifted objects are placed at their comoving distance. --- */
const deepList = [];
/* colour by lookback time: recent light is blue-white, the oldest is red */
function lookbackColor(gyr) {
  const S = [[0, [200, 222, 255]], [3, [255, 244, 214]], [7, [255, 196, 122]], [11, [255, 132, 96]], [13.8, [255, 90, 90]]];
  let i = 0; while (i < S.length - 2 && gyr > S[i + 1][0]) i++;
  const t = clamp((gyr - S[i][0]) / (S[i + 1][0] - S[i][0]), 0, 1);
  const c = S[i][1].map(function (v, k) { return Math.round(lerp(v, S[i + 1][1][k], t)); });
  return '#' + c.map(function (v) { return v.toString(16).padStart(2, '0'); }).join('');
}
function hexRgb(h) { const n = parseInt(h.slice(1), 16); return 'rgb(' + (n >> 16 & 255) + ',' + (n >> 8 & 255) + ',' + (n & 255) + ')'; }
(function () {
  if (typeof DEEP_UNIVERSE === 'undefined') return;
  const MLY_AU = 1e6 * LY_AU;
  const SUB = { quasar: 'Quasar', blazar: 'Blazar', lensed: 'Lensed quasar', cluster: 'Galaxy cluster', grb: 'Gamma-ray burst',
                gw: 'Gravitational-wave source', frb: 'Fast radio burst', early: 'Early galaxy', star: 'Most distant star' };
  const rnd = mulberry32(20260929);
  const g1 = function () { return Math.sqrt(-2 * Math.log(1 - rnd() + 1e-12)) * Math.cos(TAU * rnd()); };
  DEEP_UNIVERSE.forEach(function (d) {
    const stats = {};
    let dAU;
    if (d.z != null) {
      const c = cosmo(d.z);
      dAU = c.comoving * 1e9 * LY_AU; d.lookback = c.lookback;
      stats['Redshift'] = 'z = ' + d.z;
      stats['Light left it'] = fmtGyr(c.lookback) + ' ago';
      stats['The universe was'] = fmtGyr(c.age) + ' old';
      stats['Distance now'] = fmtGly(c.comoving);
    } else {
      dAU = d.dMly * MLY_AU; d.lookback = d.dMly / 1000;
      stats['Distance'] = fmtDist(dAU);
    }
    if (d.approx) stats['Position'] = 'Approximate';
    d.stats = Object.assign(stats, d.stats || {});
    const dir = norm(raDecToVec(d.ra, d.dec, 1));
    const rLy = d.sub === 'cluster' ? (d.rMly || 4) * 1e6 : (d.sub === 'quasar' || d.sub === 'blazar' || d.sub === 'lensed') ? 50000 : d.sub === 'early' ? 2000 : 30000;
    const o = addObject({
      id: d.id, name: d.name, type: SUB[d.sub] + (d.z != null && d.z >= 1 ? ' · z ' + d.z : ''), kind: 'distant', sub: d.sub,
      category: 'deepuniverse', color: d.color || lookbackColor(d.lookback), radiusAU: rLy * LY_AU, dir: dir,
      fixed: { x: dir.x * dAU, y: dir.y * dAU, z: dir.z * dAU }, data: d, prio: d.prio || 45, minPx: 2.5
    });
    o.pos = o.fixed;
    if (d.sub === 'cluster') {                    /* member galaxies, fixed in 3D so the cluster has depth */
      o.members = [];
      for (let k = 0; k < (d.members || 80); k++) {
        const s = k === 0 ? 0 : 0.42;
        o.members.push({ x: g1() * s, y: g1() * s, z: g1() * s, b: k === 0 ? 1 : Math.pow(rnd(), 2.5) });
      }
    }
    if (d.jet === 'us') o.jetDir = { x: -dir.x, y: -dir.y, z: -dir.z };         /* a blazar's jet points at us */
    else if (d.jet) o.jetDir = norm(cross(dir, { x: 0, y: 0, z: 1 }));        /* seen side-on */
    deepList.push(o);
  });
})();
(function () {
  if (typeof LARGE_STRUCTURES === 'undefined' || !STRUCT) return;
  const MLY_AU = 1e6 * LY_AU;
  const nrm = galacticToEcl(47.37 * DEG, 6.32 * DEG);
  const viewFrom = function (c) {
    const toUs = norm({ x: -c.x, y: -c.y, z: -c.z });
    const dd = toUs.x * nrm.x + toUs.y * nrm.y + toUs.z * nrm.z;
    const u = norm({ x: toUs.x - dd * nrm.x, y: toUs.y - dd * nrm.y, z: toUs.z - dd * nrm.z });
    const cA = Math.cos(38 * DEG), sA = Math.sin(38 * DEG);
    const v = norm({ x: u.x * cA + nrm.x * sA, y: u.y * cA + nrm.y * sA, z: u.z * cA + nrm.z * sA });
    return { yaw: Math.atan2(v.y, v.x), pitch: Math.asin(clamp(v.z, -1, 1)) };
  };
  const TYPE = { wall: 'Wall of galaxies', void: 'Cosmic void', supercluster: 'Supercluster' };
  LARGE_STRUCTURES.forEach(function (s) {
    let centre, R, path = null;
    if (s.path) {
      path = s.path.map(function (q) { return raDecToVec(q[0], q[1], s.dMly * MLY_AU); });
      centre = { x: 0, y: 0, z: 0 };
      path.forEach(function (q) { centre.x += q.x / path.length; centre.y += q.y / path.length; centre.z += q.z / path.length; });
      R = 0; path.forEach(function (q) { R = Math.max(R, Math.hypot(q.x - centre.x, q.y - centre.y, q.z - centre.z)); });
    } else {
      centre = raDecToVec(s.ra, s.dec, s.dMly * MLY_AU); R = s.rMly * MLY_AU;
    }
    s.stats = Object.assign({ 'Distance': fmtDist(Math.hypot(centre.x, centre.y, centre.z)) }, s.stats || {});
    const o = addObject({
      id: s.id, name: s.name, type: TYPE[s.kind], kind: 'region', category: 'deepuniverse', color: s.color,
      radiusAU: R, regionRadius: R, fixed: { x: centre.x, y: centre.y, z: centre.z }, data: s, prio: 50, minPx: 0, noRender: true,
      centre: true, viewFill: 0.7, view: viewFrom(centre)
    });
    o.pos = o.fixed;
    let cloud = null;
    if (path) {                                 /* a filament of galaxies scattered along the wall */
      const rnd = mulberry32(s.id.length * 977 + 13);
      const g1 = function () { return Math.sqrt(-2 * Math.log(1 - rnd() + 1e-12)) * Math.cos(TAU * rnd()); };
      const segs = []; let tot = 0;
      for (let k = 1; k < path.length; k++) { const L = Math.hypot(path[k].x - path[k - 1].x, path[k].y - path[k - 1].y, path[k].z - path[k - 1].z); segs.push(L); tot += L; }
      const n = 420, w = R * 0.07;
      cloud = { n: n, x: new Float64Array(n), y: new Float64Array(n), z: new Float64Array(n), b: new Float32Array(n) };
      for (let k = 0; k < n; k++) {
        let r = rnd() * tot, sI = 0; while (sI < segs.length - 1 && r > segs[sI]) { r -= segs[sI]; sI++; }
        const t = r / segs[sI], A = path[sI], B = path[sI + 1];
        const clump = rnd() < 0.3 ? 0.35 : 1;
        cloud.x[k] = lerp(A.x, B.x, t) + g1() * w * clump; cloud.y[k] = lerp(A.y, B.y, t) + g1() * w * clump; cloud.z[k] = lerp(A.z, B.z, t) + g1() * w * 0.5 * clump;
        cloud.b[k] = Math.pow(rnd(), 2);
      }
    }
    STRUCT.regions.push({ obj: o, centre: o.pos, R: R, normal: nrm, col: s.color, path: path, cloud: cloud, isVoid: s.kind === 'void' });
  });
})();
/* which galaxies belong to the Local Group (within 5 million light years of its centre) */
if (STRUCT) objects.forEach(function (o) {
  if (o.category !== 'galaxies' || !o.pos || o.id === 'milkyway') return;
  o.lg = Math.hypot(o.pos.x - STRUCT.lgC.x, o.pos.y - STRUCT.lgC.y, o.pos.z - STRUCT.lgC.z) < 5e6 * LY_AU;
  if (o.id === 'm31' || o.id === 'm33' || o.id === 'lmc' || o.id === 'smc') o.prio = 74;   /* the big members win label space */
});
function deepVisible(o, camFromSunLy) {
  return camFromSunLy > 3e5 || selectedId === o.id || cam.follow === o.id || focusCategory === 'deepuniverse';
}
function drawDistant(o, p) {
  const range = Math.hypot(o.pos.x - camPos.x, o.pos.y - camPos.y, o.pos.z - camPos.z);
  const rPx = (o.radiusAU / range) * focal;
  const rgb = hexRgb(o.color), soft = rgb.replace('rgb(', 'rgba(').replace(')', ',');
  ctx.save();
  if (o.sub === 'cluster') {
    if (rPx > 10 && o.members) {
      const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, rPx * 0.9);
      g.addColorStop(0, o.data.bullet ? 'rgba(180,150,255,0.10)' : soft + '0.10)'); g.addColorStop(1, soft + '0)');
      ctx.globalCompositeOperation = 'lighter';
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(p.x, p.y, rPx * 0.9, 0, TAU); ctx.fill();
      if (o.data.bullet) {                  /* the two collided gas clouds (pink) and the mass that passed through (blue) */
        const blobs = [[-0.28, 0.02, '255,110,170', 0.22], [0.12, -0.02, '255,110,170', 0.18], [-0.5, 0.05, '110,150,255', 0.2], [0.36, -0.04, '110,150,255', 0.2]];
        blobs.forEach(function (b) {
          const x = p.x + b[0] * rPx, y = p.y + b[1] * rPx, r = rPx * 0.26;
          const gg = ctx.createRadialGradient(x, y, 0, x, y, r);
          gg.addColorStop(0, 'rgba(' + b[2] + ',' + b[3] + ')'); gg.addColorStop(1, 'rgba(' + b[2] + ',0)');
          ctx.fillStyle = gg; ctx.fillRect(x - r, y - r, r * 2, r * 2);
        });
      }
      const R = o.radiusAU;
      for (let k = 0; k < o.members.length; k++) {
        const m = o.members[k];
        const q = project({ x: o.pos.x + m.x * R, y: o.pos.y + m.y * R, z: o.pos.z + m.z * R });
        if (!q) continue;
        const s = k === 0 ? Math.max(2.5, rPx * 0.03) : 0.8 + m.b * 1.8;
        ctx.fillStyle = k === 0 ? 'rgba(255,236,200,0.95)' : 'rgba(255,236,205,' + (0.35 + m.b * 0.55).toFixed(2) + ')';
        ctx.beginPath(); ctx.arc(q.x, q.y, s, 0, TAU); ctx.fill();
      }
    } else {
      glow(p, 11, rgb, 0.5);
      ctx.fillStyle = '#fff4e2'; ctx.beginPath(); ctx.arc(p.x, p.y, 1.8, 0, TAU); ctx.fill();
    }
  } else if (o.sub === 'grb' || o.sub === 'gw' || o.sub === 'frb') {
    const col = o.sub === 'gw' ? '120,220,255' : o.sub === 'frb' ? '255,130,230' : '255,214,120';
    ctx.globalCompositeOperation = 'lighter';
    for (let k = 0; k < 4; k++) {
      const a = k * Math.PI / 2 + Math.PI / 4, L = 13;
      const lg = ctx.createLinearGradient(p.x, p.y, p.x + Math.cos(a) * L, p.y + Math.sin(a) * L);
      lg.addColorStop(0, 'rgba(' + col + ',0.9)'); lg.addColorStop(1, 'rgba(' + col + ',0)');
      ctx.strokeStyle = lg; ctx.lineWidth = 1.3;
      ctx.beginPath(); ctx.moveTo(p.x, p.y); ctx.lineTo(p.x + Math.cos(a) * L, p.y + Math.sin(a) * L); ctx.stroke();
    }
    ctx.strokeStyle = 'rgba(' + col + ',0.55)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(p.x, p.y, 6, 0, TAU); ctx.stroke();
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(p.x, p.y, 1.6, 0, TAU); ctx.fill();
  } else {
    ctx.globalCompositeOperation = 'lighter';
    const host = o.sub === 'quasar' || o.sub === 'blazar' || o.sub === 'lensed';
    if (host && rPx > 3) {                    /* the galaxy around the black hole */
      const g = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, rPx);
      g.addColorStop(0, 'rgba(255,240,220,0.35)'); g.addColorStop(0.35, 'rgba(220,210,235,0.12)'); g.addColorStop(1, 'rgba(200,200,240,0)');
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(p.x, p.y, rPx, 0, TAU); ctx.fill();
    }
    if (o.jetDir && rPx > 2) {                /* jets, 200,000 light years each way */
      const L = 200000 * LY_AU;
      [1, -1].forEach(function (sgn) {
        const q = project({ x: o.pos.x + o.jetDir.x * L * sgn, y: o.pos.y + o.jetDir.y * L * sgn, z: o.pos.z + o.jetDir.z * L * sgn });
        if (!q) return;
        const lg = ctx.createLinearGradient(p.x, p.y, q.x, q.y);
        lg.addColorStop(0, 'rgba(170,200,255,0.75)'); lg.addColorStop(1, 'rgba(170,200,255,0)');
        ctx.strokeStyle = lg; ctx.lineWidth = Math.max(1.5, rPx * 0.05);
        ctx.beginPath(); ctx.moveTo(p.x, p.y); ctx.lineTo(q.x, q.y); ctx.stroke();
      });
    }
    const small = o.sub === 'early' || o.sub === 'star';
    glow(p, Math.max(small ? 10 : 15, rPx * 0.25), rgb, small ? 0.55 : 0.75);
    ctx.globalCompositeOperation = 'source-over';
    /* the core takes the lookback colour: blue-white is recent light, red is the oldest */
    ctx.fillStyle = o.color; ctx.beginPath(); ctx.arc(p.x, p.y, small ? 2.4 : 3, 0, TAU); ctx.fill();
    ctx.fillStyle = 'rgba(255,255,255,0.9)'; ctx.beginPath(); ctx.arc(p.x, p.y, small ? 0.9 : 1.2, 0, TAU); ctx.fill();
  }
  ctx.restore();
}

function drawConstellationNames() {
  if (!opts.constellations || (deep.on && deep.Y > 5e4)) return;   /* the figures distort as the stars move */
  const fade = skyFade();
  if (fade <= 0.03) return;
  ctx.save();
  ctx.font = '600 10.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  ctx.textBaseline = 'middle';
  ctx.textAlign = 'center';
  if ('letterSpacing' in ctx) ctx.letterSpacing = '2.5px';
  let placed = 0;
  for (let i = 0; i < constellationList.length && placed < 8; i++) {
    const c = constellationList[i];
    const p = projectDir(c.dir);
    if (!p) continue;
    if (p.x < 70 || p.x > W - 70 || p.y < 26 || p.y > H - 26) continue;
    const w = ctx.measureText(c.label).width;
    const box = { x: p.x - w / 2 - 8, y: p.y - 9, w: w + 16, h: 18 };
    /* object labels always win — constellation names fill the empty sky */
    if (overlapsLabel(box)) continue;
    labelBoxes.push(box);
    const sel = selectedId === c.id;
    ctx.fillStyle = sel ? 'rgba(190,204,255,' + (0.95 * fade).toFixed(3) + ')'
                        : 'rgba(150,168,212,' + (0.42 * fade).toFixed(3) + ')';
    ctx.fillText(c.label, p.x, p.y);
    placed++;
  }
  ctx.restore();
}

/* ============================================================
   7. Interaction — orbit, zoom, pick
   ============================================================ */
let dragging = false, lastX = 0, lastY = 0, movedPx = 0, idle = 0;

canvas.addEventListener('pointerdown', function (e) {
  dragging = true; movedPx = 0; idle = 0;
  lastX = e.clientX; lastY = e.clientY;
  canvas.classList.add('grabbing');
  canvas.setPointerCapture(e.pointerId);
});
canvas.addEventListener('pointermove', function (e) {
  if (!dragging) return;
  const dx = e.clientX - lastX, dy = e.clientY - lastY;
  movedPx += Math.abs(dx) + Math.abs(dy);
  lastX = e.clientX; lastY = e.clientY;
  const s = 0.005 * (cam.fov / FOV_DEFAULT);       /* a magnified sky pans slower */
  if (opts.freeRotate && cam.fovGoal >= FOV_DEFAULT * 0.999) {
    cam.pendYaw -= dx * s; cam.pendPitch += dy * s;       /* any direction, over the poles */
  } else {
    cam.yawGoal -= dx * s;
    cam.pitchGoal = clamp(cam.pitchGoal + dy * s, -1.52, 1.52);   /* level horizon: stop just short of the poles */
  }
});
canvas.addEventListener('pointerup', function (e) {
  dragging = false;
  canvas.classList.remove('grabbing');
  if (movedPx < 5) pick(e.clientX, e.clientY, false);
});
canvas.addEventListener('dblclick', function (e) { pick(e.clientX, e.clientY, true); });

canvas.addEventListener('wheel', function (e) {
  e.preventDefault();
  idle = 0;
  /* deltaMode 1 = lines, 2 = pages; Firefox reports lines, so raw deltaY is tiny */
  const dy = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaMode === 2 ? e.deltaY * 400 : e.deltaY;
  zoomBy(Math.exp(clamp(dy, -120, 120) * 0.0016));
}, { passive: false });
/* Zoom means two things: moving the camera in and out, or — once you are looking at
   a telescope's patch of sky — magnifying the sky itself. Zooming out of a magnified
   view returns to the normal field of view before the camera starts moving. */
function zoomBy(f) {
  const sel = selectedId && byId[selectedId];
  const magnified = cam.fovGoal < FOV_DEFAULT * 0.999;
  if (magnified || (f < 1 && sel && sel.kind === 'deepfield')) {
    cam.fovGoal = clamp(cam.fovGoal * f, FOV_MIN, FOV_DEFAULT);
    return;
  }
  cam.distGoal = clamp(cam.distGoal * f, MIN_DIST, MAX_DIST);
}

function pick(clientX, clientY, fly) {
  const rect = canvas.getBoundingClientRect();
  const mx = clientX - rect.left, my = clientY - rect.top;
  let best = null, bestD = 26;
  for (let i = 0; i < objects.length; i++) {
    const o = objects[i];
    if (!o.screen || (o.noRender && o.kind !== 'deepfield')) continue;
    const d = Math.hypot(o.screen.x - mx, o.screen.y - my);
    const r = Math.max(o.minPx || 2, (o.radiusAU / o.screen.z) * focal);
    const reach = Math.max(14, r + 8);
    if (d < reach && d < bestD + r) { best = o; bestD = d; }
  }
  if (best) { select(best.id); if (fly) flyTo(best.id); }
  else if (!fly) { select(null); }
}

/* Distance at which an object of this radius fills `fill` of the viewport
   height. Arriving should frame the thing with space around it — you want
   to SEE Pluto, not stand on it. */
function frameDistance(radiusAU, fill) {
  return radiusAU / Math.max(Math.tan(FOV_DEFAULT / 2) * fill, 1e-12);
}

function objectViewDistance(o) {
  if (o.viewFill)             return frameDistance(o.radiusAU, o.viewFill);
  if (o.kind === 'region')    return o.regionRadius * 2.6;
  if (o.kind === 'probe') {
    if (o.mode === 'gp') return 0.00004;
    if (o.mode === 'l2' || o.mode === 'l2halo') return 0.035;          /* Earth and L2 both in view */
    const host = o.parent && byId[o.parent];
    if (host) return Math.max(frameDistance(host.radiusAU, 0.22), o.aAU ? o.aAU * 1.6 : 0);
    return 2.2;
  }
  /* a black hole's visible extent is its disc, ~9.5x the shadow radius */
  if (o.kind === 'blackhole') return frameDistance(o.radiusAU * 9.5, 0.55);
  if (o.kind === 'distant') return o.sub === 'cluster' ? frameDistance(o.radiusAU, 0.55) : frameDistance(o.jetDir ? 200000 * LY_AU : o.radiusAU * 1.4, 0.6);
  if (o.kind === 'deepsky')   return frameDistance(o.radiusAU, 0.42);
  if (o.kind === 'neutron')   return frameDistance(o.radiusAU, 0.30);
  if (o.kind === 'star')      return Math.max(frameDistance(o.radiusAU, 0.014), 0.05);
  if (o.kind === 'sun')       return frameDistance(o.radiusAU, 0.30);
  return Math.max(frameDistance(o.radiusAU, 0.32), 1.2e-5);
}

/* Sun-lit approach. Bodies shine only by reflected sunlight, so arriving on
   the far side means staring at an unlit disc against black sky. Sit on the
   sunward side, swung ~30 deg off-axis so a terminator still shows. */
function sunlitAngles(o) {
  const L = (o.lightId && byId[o.lightId]) ? byId[o.lightId].pos : { x: 0, y: 0, z: 0 };
  const dx = L.x - o.pos.x, dy = L.y - o.pos.y, dz = L.z - o.pos.z;
  const m = Math.hypot(dx, dy, dz);
  if (!isFinite(m) || m < 1e-9) return null;
  const s = { x: dx / m, y: dy / m, z: dz / m };   /* body -> its star */
  let perp = cross(s, { x: 0, y: 0, z: 1 });
  if (Math.hypot(perp.x, perp.y, perp.z) < 1e-6) perp = { x: 1, y: 0, z: 0 };
  perp = norm(perp);
  const k = Math.tan(30 * DEG);
  const u = norm({ x: s.x + perp.x * k, y: s.y + perp.y * k, z: s.z + perp.z * k + 0.30 });
  return { yaw: Math.atan2(u.y, u.x), pitch: clamp(Math.asin(clamp(u.z, -1, 1)), -1.2, 1.2) };
}

function flyTo(id) {
  const o = byId[id];
  if (!o) return;
  if (o.kind === 'constellation' || o.kind === 'deepfield' || o.kind === 'ngc') {
    /* a direction, not a place — turn to look at it, and magnify a patch so it fills the view */
    const ang = o.angArcmin || 0;
    lookAtDir(o.dir, ang > 0 ? (ang / 0.45) / 60 * DEG : FOV_DEFAULT);
    pushHistory(id);
    toast(o.kind === 'deepfield' ? 'Looking at ' + o.name + ' · scroll to zoom the view' : 'Looking toward ' + o.name);
    return;
  }
  cam.fovGoal = FOV_DEFAULT;       /* flying somewhere always leaves the telescope view */
  if (o.kind === 'region') {
    cam.follow = null;
    /* belts are centred on the Sun; a galaxy has its own centre and its own best angle */
    const dest = o.centre ? o : { pos: { x: 0, y: 0, z: 0 }, id: 'sun' };
    let yaw = cam.yawGoal, pitch = 0.62;
    if (o.view) { let dy = o.view.yaw - cam.yawGoal; dy -= TAU * Math.round(dy / TAU); yaw = cam.yawGoal + dy; pitch = o.view.pitch; }
    cam.flight = {
      t: 0, dur: o.centre ? 3.4 : 2.0, obj: dest,
      fx: cam.tx, fy: cam.ty, fz: cam.tz,
      fd: cam.dist, td: objectViewDistance(o),
      fyaw: cam.yawGoal, fpitch: cam.pitchGoal, yaw: yaw, pitch: pitch,
      noFollow: true          /* a region is a volume, not a body to lock onto */
    };
    pushHistory(id);
    toast('Showing ' + o.name);
    return;
  }
  const d0 = Math.hypot(o.pos.x - cam.tx, o.pos.y - cam.ty, o.pos.z - cam.tz);
  const td = objectViewDistance(o);
  const dur = clamp(1.3 + Math.log10(1 + d0 / 5) * 0.5 + Math.abs(Math.log10(td / cam.dist)) * 0.28, 1.3, 3.6);
  /* sunlit approach for anything that shines by reflected light */
  const lit = (o.kind === 'planet' || o.kind === 'dwarf' || o.kind === 'moon' ||
               o.kind === 'exoplanet') ? sunlitAngles(o) : null;
  let yawGoal = cam.yawGoal;
  let pitchGoal = o.kind === 'blackhole' ? 0.22 : 0.32;
  if (lit) { yawGoal = lit.yaw; pitchGoal = lit.pitch; }
  /* rotate the short way round rather than the long way */
  let dy = yawGoal - cam.yawGoal;
  dy -= TAU * Math.round(dy / TAU);

  cam.follow = null;
  updateFocusChip();
  cam.flight = {
    t: 0, dur: dur, obj: o,
    fx: cam.tx, fy: cam.ty, fz: cam.tz,
    fd: cam.dist, td: td,
    fyaw: cam.yawGoal, fpitch: cam.pitchGoal,
    yaw: cam.yawGoal + dy, pitch: pitchGoal
  };
  pushHistory(id);
  toast('Travelling to ' + o.name);
}

/* ============================================================
   8. Interface
   ============================================================ */
const $ = function (id) { return document.getElementById(id); };
function esc(t) {
  return String(t).replace(/[&<>"']/g, function (c) {
    return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
  });
}

const ICONS = {
  star: '<polygon points="12 2 15 9 22 9.5 16.5 14.2 18.3 21 12 17.1 5.7 21 7.5 14.2 2 9.5 9 9"/>',
  orbit: '<circle cx="12" cy="12" r="3.2"/><ellipse cx="12" cy="12" rx="10" ry="4.2"/><ellipse cx="12" cy="12" rx="10" ry="4.2" transform="rotate(60 12 12)"/>',
  moon: '<path d="M17 3a9 9 0 1 0 4 8 7 7 0 0 1-4-8z"/>',
  dwarf: '<circle cx="12" cy="12" r="6"/><circle cx="19.5" cy="6" r="1.5"/>',
  belt: '<ellipse cx="12" cy="12" rx="9.5" ry="3.6"/><circle cx="12" cy="12" r="2"/><circle cx="20" cy="13.4" r=".8" fill="currentColor"/><circle cx="4.5" cy="10.6" r=".8" fill="currentColor"/>',
  rocket: '<path d="M12 2.5c2.5 2 4 6 4 10 0 2-.5 4-1 5.5l-3 2-3-2c-.5-1.5-1-3.5-1-5.5 0-4 1.5-8 4-10z"/><path d="M9.5 15 6 16.5 7 19"/><path d="M14.5 15 18 16.5 17 19"/><circle cx="12" cy="10" r="1.6"/>',
  sparkle: '<path d="M12 3v6M12 15v6M3 12h6M15 12h6"/><circle cx="12" cy="12" r="2.2"/>',
  nebula: '<path d="M6.5 17a4 4 0 0 1-1-7.9A5 5 0 0 1 15 7.1 4.5 4.5 0 0 1 17.5 16z"/><circle cx="10" cy="12" r="1" fill="currentColor" stroke="none"/><circle cx="14" cy="13.5" r="1" fill="currentColor" stroke="none"/>',
  galaxy: '<path d="M12 4c4.4 0 8 3.6 8 8"/><path d="M12 20c-4.4 0-8-3.6-8-8"/><path d="M12 8.5a3.5 3.5 0 0 1 3.5 3.5"/><path d="M12 15.5A3.5 3.5 0 0 1 8.5 12"/><circle cx="12" cy="12" r="1.2" fill="currentColor" stroke="none"/>',
  telescope: '<path d="M3.5 14.5 14 8.5l3 5.2-10.5 6z"/><path d="M14 8.5 16.4 4l3.6 2.1-2.9 4.6"/><path d="M9 17.5 7.5 21"/><path d="M13 15.5 15 21"/>',
  quasar: '<circle cx="12" cy="12" r="2.6" fill="currentColor" stroke="none"/><circle cx="12" cy="12" r="6.5" stroke-dasharray="2 2.5"/><path d="M12 1.5v5M12 17.5v5"/>',
  address: '<circle cx="12" cy="12" r="9"/><circle cx="12" cy="12" r="4.6"/><circle cx="12" cy="12" r="1.3" fill="currentColor" stroke="none"/>',
  deepfield: '<rect x="4" y="4" width="16" height="16" rx="1" stroke-dasharray="4 3"/><circle cx="9.5" cy="10" r="1" fill="currentColor" stroke="none"/><circle cx="14" cy="13.5" r="1.3" fill="currentColor" stroke="none"/><circle cx="15" cy="8.5" r="0.8" fill="currentColor" stroke="none"/>',
  constellation: '<circle cx="6" cy="8" r="1.2" fill="currentColor" stroke="none"/><circle cx="12" cy="5.5" r="1" fill="currentColor" stroke="none"/><circle cx="17" cy="10" r="1.3" fill="currentColor" stroke="none"/><circle cx="9" cy="15" r="1.1" fill="currentColor" stroke="none"/><circle cx="15.5" cy="18" r="1" fill="currentColor" stroke="none"/>',
  satellite: '<rect x="9.5" y="9.5" width="5" height="5" rx="1"/><path d="M3 8l4 4-2 2-4-4zM21 16l-4-4 2-2 4 4zM9.5 12H7M17 12h-2.5"/>',
  asteroid: '<path d="M8 3.5 5 6l-1.5 4L5 15l3.5 3 4.5 1.5 4-2 2.5-3.5-1-4.5-3-4L10 4z"/><circle cx="10.5" cy="9.5" r="1" fill="currentColor" stroke="none"/><circle cx="14" cy="13" r="1.2" fill="currentColor" stroke="none"/>',
  comet: '<circle cx="7" cy="17" r="3"/><path d="M9.5 14.5 20 4M10.5 17 19 12M7 13.5 12 5"/>',
  exoplanet: '<circle cx="12" cy="12" r="4.6"/><ellipse cx="12" cy="12" rx="10" ry="3.6" transform="rotate(-25 12 12)"/>',
  remnant: '<circle cx="12" cy="12" r="2.6" fill="currentColor" stroke="none"/><path d="M12 2.5v4M12 17.5v4M4.4 7.2l3.4 2M16.2 14.8l3.4 2M4.4 16.8l3.4-2M16.2 9.2l3.4-2"/>',
  cluster: '<circle cx="9" cy="9" r="1.3" fill="currentColor" stroke="none"/><circle cx="15" cy="8" r="1" fill="currentColor" stroke="none"/><circle cx="12" cy="12.5" r="1.5" fill="currentColor" stroke="none"/><circle cx="8" cy="15" r="1.1" fill="currentColor" stroke="none"/><circle cx="16" cy="15.5" r="1.2" fill="currentColor" stroke="none"/>',
  blackhole: '<ellipse cx="12" cy="12" rx="10" ry="3.4"/><circle cx="12" cy="12" r="4" fill="currentColor" stroke="none"/>'
};
function svg(name) {
  return '<svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" ' +
         'stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">' + (ICONS[name] || '') + '</svg>';
}

/* ---- browse panel ---- */
const panel = $('panel'), panelRoot = $('panelRoot'), panelDetail = $('panelDetail');

CATEGORIES.forEach(function (c) {
  const li = document.createElement('li');
  li.innerHTML =
    '<button class="category-btn" data-cat="' + c.id + '">' +
      '<span class="category-icon cat-' + c.id + '">' + svg(c.icon) + '</span>' +
      '<span class="category-text"><strong>' + c.name + '</strong><span>' + c.sub + '</span></span>' +
      '<span class="category-chev"><svg viewBox="0 0 24 24" width="15" height="15" fill="none" ' +
        'stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">' +
        '<polyline points="9 6 15 12 9 18"/></svg></span>' +
    '</button>';
  li.querySelector('button').addEventListener('click', function () { openCategory(c.id); });
  $('categoryList').appendChild(li);
});

function categoryMembers(catId) {
  if (catId === 'featured') return FEATURED_IDS.map(function (id) { return byId[id]; }).filter(Boolean);
  const list = objects.filter(function (o) { return o.category === catId; });
  if (catId === 'galaxies' && byId.milkyway) return [byId.milkyway].concat(list.filter(function (o) { return o.id !== 'milkyway'; }));
  if (catId === 'structure') return ['arm-orion', 'milkyway', 'localgroup', 'virgosc', 'laniakea', 'observable', 'bar', 'arm-sagittarius', 'arm-perseus', 'arm-scutum', 'arm-norma']
    .map(function (id) { return byId[id]; }).filter(Boolean);
  return list;
}

function openCategory(catId) {
  const cat = CATEGORIES.find(function (c) { return c.id === catId; });
  const list = categoryMembers(catId);
  $('detailTitle').textContent = cat.name;
  $('detailSub').textContent = cat.sub;
  const ul = $('detailList');
  ul.innerHTML = '';
  list.forEach(function (o) {
    const li = document.createElement('li');
    const badge = o.kind === 'deepfield' ? (o.data.planned ? 'planned' : String(o.data.year))
                : o.kind === 'constellation' ? o.data.area.toLocaleString() + ' sq°'
                : o.kind === 'star' ? (o.data.d < 100 ? o.data.d.toFixed(2) + ' ly' : Math.round(o.data.d).toLocaleString() + ' ly')
                : o.kind === 'blackhole' ? Math.round(o.data.d).toLocaleString() + ' ly'
                : o.kind === 'distant' ? (o.data.z != null ? 'z ' + o.data.z : o.data.dMly + ' Mly') : '';
    li.innerHTML =
      '<button class="object-btn">' +
        '<span class="object-dot" style="background:' + o.color + ';color:' + o.color + '"></span>' +
        '<span class="object-text"><strong>' + o.name + '</strong><span>' + o.type + '</span></span>' +
        (badge ? '<span class="object-badge">' + badge + '</span>' : '') +
      '</button>';
    li.querySelector('button').addEventListener('click', function () { select(o.id); flyTo(o.id); });
    ul.appendChild(li);
  });
  focusCategory = catId;
  panelRoot.classList.add('hidden');
  panelDetail.classList.remove('hidden');
}
$('detailBack').addEventListener('click', function () {
  focusCategory = null;
  panelDetail.classList.add('hidden');
  panelRoot.classList.remove('hidden');
});
$('panelClose').addEventListener('click', function () {
  panel.classList.add('hidden');
  $('reopenBtn').classList.add('show');
});
$('reopenBtn').addEventListener('click', function () {
  panel.classList.remove('hidden');
  $('reopenBtn').classList.remove('show');
});


/* ---- Related objects: what orbits, lands on, or circles the selected body ----
   Shown in the info panel as clickable chips. Earth additionally offers the
   artificial-satellite layer as an opt-in, so 16,000 dots never appear uninvited. */
function probeHost(o) {
  if (o.kind !== 'probe') return null;
  if (o.parent) return o.parent;
  if (o.mode === 'l2' || o.mode === 'l2halo' || o.mode === 'gp') return 'earth';
  return null;
}
function relatedGroups(o) {
  const groups = [];
  const push = function (title, list, badge) { if (list.length) groups.push({ title: title, list: list, badge: badge }); };
  const yearBadge = function (x) { return x.data.planned ? 'planned' : String(x.data.year); };
  if (o.category === 'telescopes' && o.data && o.data.scope) {
    const seen = deepFieldList.filter(function (d) { return d.data.scopes.indexOf(o.data.scope) >= 0; })
      .sort(function (p, q) { return p.data.year - q.data.year; });
    push('What it saw', seen, yearBadge);
  }
  if (o.id === 'localgroup') {
    const mem = [byId.milkyway].concat(objects.filter(function (x) { return x.lg; })
      .sort(function (p, q) { return (q.radiusAU - p.radiusAU); })).filter(Boolean);
    push('Member galaxies', mem);
  }
  if (o.id === 'milkyway' || o.id === 'm31') push('Satellite galaxies', objects.filter(function (x) { return x.data && x.data.within === o.id && x.category === 'galaxies'; }));
  if (o.kind === 'deepfield') {
    push('Seen by', objects.filter(function (x) { return x.category === 'telescopes' && x.data && o.data.scopes.indexOf(x.data.scope) >= 0; }));
    push('Same patch of sky', deepFieldList.filter(function (d) {
      return d !== o && angSepArcmin(d.dir, o.dir) < Math.max(20, 0.6 * (d.data.arcmin + o.data.arcmin));
    }), yearBadge);
  }
  if (o.data && (o.data.links || o.data.link)) {
    const links = (o.data.links || [o.data.link]).map(function (id) { return byId[id]; }).filter(Boolean);
    push('In the map', links);
  }
  if (o.kind === 'sun') {
    push('Planets', objects.filter(function (x) { return x.kind === 'planet'; }));
    push('Dwarf planets', objects.filter(function (x) { return x.kind === 'dwarf' && x.category === 'dwarfs'; }));
  }
  push('Moons', objects.filter(function (x) { return x.kind === 'moon' && x.parent === o.id; }));
  push('Planets', objects.filter(function (x) { return x.kind === 'exoplanet' && x.parent === o.id; }));
  push('Spacecraft', objects.filter(function (x) { return x.kind === 'probe' && x.category !== 'satellites' && probeHost(x) === o.id; }));
  const host = (o.kind === 'moon' || o.kind === 'exoplanet') ? o.parent : probeHost(o);
  if (host && byId[host]) push('Orbits', [byId[host]]);
  return groups;
}
function chipHtml(x, badge) {
  return '<button class="chip" data-id="' + esc(x.id) + '" style="--c:' + esc(x.color || '#aaa') + '">' + esc(x.name) +
         (badge ? '<span class="chip-badge">' + esc(badge) + '</span>' : '') + '</button>';
}
function renderRelated(o) {
  const box = $('infoRelated'); if (!box) return;
  let html = '';
  relatedGroups(o).forEach(function (g) {
    html += '<div class="rel-group"><div class="rel-title">' + esc(g.title) +
            (g.list.length > 1 ? ' <span>' + g.list.length + '</span>' : '') + '</div><div class="rel-chips">' +
            g.list.slice(0, 40).map(function (x) { return chipHtml(x, g.badge ? g.badge(x) : null); }).join('') + '</div>' +
            (g.title === 'What it saw' ? '<div class="rel-note">Pick a view to turn toward that patch of sky and zoom in on it.</div>' : '') + '</div>';
  });
  if (o.kind === 'deepfield' && o.data.highlights && o.data.highlights.length) {
    html += '<div class="rel-group"><div class="rel-title">Notable finds</div>';
    o.data.highlights.forEach(function (h, i) {
      html += '<div class="find"><div class="find-head"><span class="find-name">' + esc(h.name) + '</span>' +
              (h.z != null ? '<span class="find-z">z = ' + h.z + '</span>' : '') +
              (h.dir ? '<button class="btn tiny find-show" data-i="' + i + '">Show</button>' : '') + '</div>' +
              (h.z != null ? '<div class="find-meta">' + esc(cosmoSentence(h.z)) + ' Its light is stretched ' + (1 + h.z).toFixed(1) + '× by that expansion.</div>' : '') +
              (h.note ? '<div class="find-note">' + esc(h.note) + '</div>' : '') + '</div>';
    });
    html += '</div>';
  }
  const S = (typeof bulkSat !== 'undefined') ? bulkSat : null;
  if (o.id === 'earth' && S) {
    const feat = objects.filter(function (x) { return x.category === 'satellites'; });
    html += '<div class="rel-group"><div class="rel-title">Artificial satellites <span>' + S.n.toLocaleString() + ' tracked</span></div>' +
            '<button class="btn tiny rel-toggle' + (opts.satellites ? ' on' : '') + '" id="btnSatLayer">' +
            (opts.satellites ? 'Hide the satellite layer' : 'Show all ' + S.n.toLocaleString() + ' satellites') + '</button>' +
            '<div class="rel-chips">' + feat.map(chipHtml).join('') + '</div></div>';
  }
  const addr = cosmicAddress(o);
  if (addr.length) {
    html += '<div class="rel-group"><div class="rel-title">Cosmic address</div><div class="address"><span class="addr-here">' + esc(o.name) + '</span>' +
      addr.map(function (s) {
        return '<span class="addr-sep">›</span>' + (s.id && byId[s.id]
          ? '<button class="chip addr" data-id="' + esc(s.id) + '" style="--c:' + esc(byId[s.id].color || '#aaa') + '">' + esc(s.name) + '</button>'
          : '<span class="addr-plain">' + esc(s.name) + '</span>');
      }).join('') + '</div></div>';
  }
  box.innerHTML = html;
  box.classList.toggle('hidden', !html);
  box.querySelectorAll('.chip').forEach(function (b) {
    b.addEventListener('click', function () { const id = b.getAttribute('data-id'); select(id); flyTo(id); });
  });
  box.querySelectorAll('.find-show').forEach(function (b) {
    b.addEventListener('click', function () {
      const i = +b.getAttribute('data-i'), h = o.data.highlights[i];
      if (!h || !h.dir) return;
      o.data._shown = i;
      lookAtDir(h.dir, 2.5 / 60 * DEG);
      toast('Centred on ' + h.name);
    });
  });
  const t = $('btnSatLayer');
  if (t) t.addEventListener('click', function () {
    opts.satellites = !opts.satellites;
    const cb = $('optSatellites'); if (cb) cb.checked = opts.satellites;
    renderRelated(o);
  });
}

/* ---- info panel ---- */
function select(id) {
  selectedId = id;
  const info = $('infoPanel');
  if (!id || !byId[id]) { info.classList.add('hidden'); return; }
  const o = byId[id];
  info.classList.remove('hidden');
  $('infoName').textContent = o.name;
  $('infoType').textContent = o.type;
  const sw = $('infoSwatch');
  sw.style.background = o.kind === 'blackhole' ? '#000' : o.color;
  sw.style.color = o.color;
  $('infoDesc').textContent = (o.data && o.data.desc) || '';
  renderStats(o);
  renderRelated(o);
  $('btnFollow').classList.toggle('on', cam.follow === id);
  $('btnCompare').style.display = sizable(o) ? '' : 'none';
  refreshGotoLabel();
}

function renderStats(o) {
  const dl = $('infoStats');
  dl.innerHTML = '';
  const rows = [];
  const stats = (o.data && o.data.stats) || {};
  Object.keys(stats).forEach(function (k) { rows.push([k, stats[k]]); });

  if (o.kind !== 'region' && o.kind !== 'constellation' && o.kind !== 'deepfield' && o.kind !== 'ngc' && o.kind !== 'distant') {
    const rSun = Math.hypot(o.pos.x, o.pos.y, o.pos.z);
    if (o.kind === 'star' || o.kind === 'blackhole') {
      rows.push(['Distance from Sun', fmtDist(rSun)]);
    } else if (o.kind !== 'sun') {
      rows.push(['Distance from Sun', fmtDist(rSun)]);
      const e = byId.earth;
      rows.push(['Distance from Earth', fmtDist(Math.hypot(o.pos.x - e.pos.x, o.pos.y - e.pos.y, o.pos.z - e.pos.z))]);
    }
    if (o.kind === 'star' && o.absMag != null) {
      rows.push(['Apparent magnitude', o.appMag.toFixed(2)]);
      rows.push(['Absolute magnitude', o.absMag.toFixed(2)]);
    }
    const lightMin = Math.hypot(o.pos.x, o.pos.y, o.pos.z) * 8.317;
    if (o.kind !== 'sun' && lightMin < 60 * 24) {
      rows.push(['Light from the Sun', lightMin < 60 ? lightMin.toFixed(1) + ' min' : (lightMin / 60).toFixed(1) + ' hours']);
    }
  }
  rows.forEach(function (r) {
    const dt = document.createElement('dt'); dt.textContent = r[0];
    const dd = document.createElement('dd'); dd.textContent = r[1];
    dl.appendChild(dt); dl.appendChild(dd);
  });
}

$('infoClose').addEventListener('click', function () { select(null); });
$('btnGoto').addEventListener('click', function () { if (selectedId) flyTo(selectedId); });
/* the action reads differently for a direction than for a destination */
function refreshGotoLabel() {
  const o = selectedId && byId[selectedId];
  const aimOnly = o && (o.kind === 'constellation' || o.kind === 'deepfield' || o.kind === 'ngc');
  $('btnGoto').textContent = aimOnly ? 'Look here' : 'Fly here';
  $('btnFollow').style.display = aimOnly ? 'none' : '';
}
$('btnFollow').addEventListener('click', function () {
  if (!selectedId) return;
  cam.follow = (cam.follow === selectedId) ? null : selectedId;
  $('btnFollow').classList.toggle('on', cam.follow === selectedId);
  updateFocusChip();
});

function updateFocusChip() {
  const chip = $('focusChip');
  if (cam.follow && byId[cam.follow]) {
    chip.classList.remove('hidden');
    $('focusName').textContent = 'Following ' + byId[cam.follow].name;
  } else {
    chip.classList.add('hidden');
  }
  if (selectedId) $('btnFollow').classList.toggle('on', cam.follow === selectedId);
}
$('focusClear').addEventListener('click', function () { cam.follow = null; updateFocusChip(); });

/* Catalogue rows are packed arrays, not objects. When one is searched for or
   selected, promote just that row into a normal object so the rest of the app
   (info panel, fly-to, history) needs no special cases. */
function bulkStarObject(i) {
  const id = 'hyg' + i;
  if (byId[id]) return byId[id];
  const S = bulkStars, e = S.named[i];
  const dAU = Math.hypot(S.x[i], S.y[i], S.z[i]);
  const ly = dAU / LY_AU, mag = S.mag[i] * 0.01;
  const spect = e ? e[3] : '';
  const sp = SPEC[(spect || 'G').charAt(0)] || SPEC.G;
  const o = addObject({
    id: id, name: e ? e[1] : 'HYG ' + i, type: 'Star' + (spect ? ' · ' + spect : ''),
    kind: 'star', category: 'stars', color: ciColor(S.ci[i] * 0.001),
    radiusAU: sp.r * R_SUN_AU, appMag: mag,
    absMag: mag - 5 * (Math.log(ly / (LY_AU && PC_LY) || 1) / Math.LN10),
    fixed: { x: S.x[i], y: S.y[i], z: S.z[i] }, prio: 38 - mag, minPx: 1.2,
    data: { d: ly, desc: 'From the HYG catalogue of nearby and naked-eye stars.',
            stats: { 'Magnitude': mag.toFixed(2), 'Distance': fmtLy(ly),
                     'Spectral type': spect || '—', 'Constellation': e ? e[2] : '—' } }
  });
  o.absMag = mag - 5 * (Math.log(ly / PC_LY) / Math.LN10) + 5;
  o.pos = o.fixed;
  return o;
}

function bulkDsoObject(i) {
  const id = 'ngc_' + i;
  if (byId[id]) return byId[id];
  const D = bulkDso, INV = 1 / 32767;
  const dir = norm({ x: D.x[i] * INV, y: D.y[i] * INV, z: D.z[i] * INV });
  const mag = D.mag[i] * 0.01, maj = D.maj[i] * 0.01;
  const tn = D.typeNames[D.type[i]];
  const full = { G: 'Galaxy', GPair: 'Galaxy pair', GTrpl: 'Galaxy triplet', GGroup: 'Galaxy group',
                 PN: 'Planetary nebula', OCl: 'Open cluster', GCl: 'Globular cluster',
                 'Cl+N': 'Cluster with nebulosity', HII: 'H II region', Neb: 'Nebula',
                 RfN: 'Reflection nebula', SNR: 'Supernova remnant', EmN: 'Emission nebula',
                 DrkN: 'Dark nebula', '*Ass': 'Stellar association' }[tn] || tn;
  const o = addObject({
    id: id, name: D.common[i] || D.names[i], type: full, kind: 'ngc', category: 'deepsky',
    color: '#cdd6ee', radiusAU: 0, dir: dir,
    fixed: { x: dir.x * 1e11, y: dir.y * 1e11, z: dir.z * 1e11 },
    prio: 25, minPx: 0, noRender: true, angArcmin: maj > 0 ? maj : 0,
    data: { desc: 'From the OpenNGC catalogue. No reliable distance is published for most NGC entries, so this is shown as a direction on the sky rather than placed in depth.',
            stats: { 'Catalogue name': D.names[i], 'Type': full,
                     'Magnitude': mag > 90 ? '—' : mag.toFixed(1),
                     'Apparent size': maj > 0 ? maj.toFixed(1) + ' arcmin' : '—' } }
  });
  o.pos = o.fixed;
  return o;
}

function bulkExoObject(i) {
  const id = 'exo_' + i;
  if (byId[id]) return byId[id];
  const E = bulkExo, INV = 1 / 32767;
  const pc = E.dist[i] / 1000, dAU = pc * PC_AU;
  const dir = norm({ x: E.x[i] * INV, y: E.y[i] * INV, z: E.z[i] * INV });
  const rE = E.rad[i] * 0.01, a = E.a[i] / 100000, per = E.per[i] * 0.01;
  const o = addObject({
    id: id, name: E.names[i], type: 'Exoplanet · ' + E.hosts[i], kind: 'exoplanet',
    category: 'exoplanets', color: '#9fb4c8', color2: '#9fb4c8',
    radiusAU: (Math.max(rE, 0.4) * 6371) / AU_KM,
    fixed: { x: dir.x * dAU, y: dir.y * dAU, z: dir.z * dAU }, prio: 26, minPx: 1.6,
    data: { desc: 'From the NASA Exoplanet Archive. Shown at its host system\'s position; the orbit itself is not modelled for bulk catalogue entries.',
            stats: { 'Host star': E.hosts[i], 'Distance': fmtLy(pc * PC_LY),
                     'Orbit radius': a > 0 ? a.toFixed(4) + ' AU' : '—',
                     'Year': per > 0 ? per.toFixed(2) + ' days' : '—',
                     'Radius': rE > 0 ? rE.toFixed(2) + ' × Earth' : '—',
                     'Discovered': E.year[i] || '—' } }
  });
  o.pos = o.fixed;
  return o;
}

/* Search returns lightweight rows. Promoting a catalogue entry into the live scene
   on every keystroke leaked objects into the per-frame loops and produced duplicate
   results once an entry existed in both `objects` and the catalogue. */
function bulkRow(id, name, type, color, make) {
  return { id: id, name: name, type: type, color: color, promote: make, bulkRow: true };
}
function searchBulk(q, out, limit) {
  const have = {};
  for (let i = 0; i < out.length; i++) have[out[i].id] = 1;
  if (bulkStars) {
    const nm = bulkStars.named;
    for (const k in nm) {
      if (out.length >= limit) return;
      if (nm[k][1].toLowerCase().indexOf(q) >= 0 && !have['hyg' + k])
        out.push(bulkRow('hyg' + k, nm[k][1], 'Star', '#dfe6ff', bulkStarObject.bind(null, +k)));
    }
  }
  if (bulkDso) {
    for (let i = 0; i < bulkDso.n && out.length < limit; i++) {
      if ((bulkDso.lname[i].indexOf(q) >= 0 || (bulkDso.lcommon[i] && bulkDso.lcommon[i].indexOf(q) >= 0))
          && !have['ngc_' + i])
        out.push(bulkRow('ngc_' + i, bulkDso.common[i] || bulkDso.names[i],
                         bulkDso.typeNames[bulkDso.type[i]], '#cdd6ee', bulkDsoObject.bind(null, i)));
    }
  }
  if (bulkGal) {
    for (let i = 0; i < bulkGal.n && out.length < limit; i++) {
      if (bulkGal.lname[i].indexOf(q) >= 0 && !have['gal' + i]) {
        const gt = galType(bulkGal.t[i]);
        out.push(bulkRow('gal' + i, bulkGal.names[i], gt[0], gt[1], bulkGalObject.bind(null, i)));
      }
    }
  }
  if (bulkSat) {
    for (let i = 0; i < bulkSat.n && out.length < limit; i++) {
      if (bulkSat.lname[i].indexOf(q) >= 0 && !have['sat' + i]) {
        const cls = SAT_CLASS[bulkSat.cls[i]];
        out.push(bulkRow('sat' + i, bulkSat.names[i], cls.label, cls.color, bulkSatObject.bind(null, i)));
      }
    }
  }
  if (bulkSbdb) {
    for (let i = 0; i < bulkSbdb.n && out.length < limit; i++) {
      if ((bulkSbdb.lname[i].indexOf(q) >= 0 || (bulkSbdb.lkey[i] && bulkSbdb.lkey[i].indexOf(q) >= 0)) && !have['sb' + i]) {
        const cls = SBDB_CLASS[bulkSbdb.cls[i]];
        out.push(bulkRow('sb' + i, bulkSbdb.names[i], cls.label, cls.color, bulkSbdbObject.bind(null, i)));
      }
    }
  }
  if (bulkExo) {
    for (let i = 0; i < bulkExo.n && out.length < limit; i++) {
      if (bulkExo.lname[i].indexOf(q) >= 0 && !have['exo_' + i])
        out.push(bulkRow('exo_' + i, bulkExo.names[i], 'Exoplanet · ' + bulkExo.hosts[i],
                         '#9fb4c8', bulkExoObject.bind(null, i)));
    }
  }
}


/* ---- Trip planner: pick two places, see the live distance and how long the
   journey takes at real vehicle speeds, plus the minimum-energy transfer. ---- */
const route = { from: 'earth', to: null, open: false };
const MU_SUN = 2.9591220828e-4;                  /* AU^3 / day^2 */
/* anything with a place in space can be a destination: planets, stars, galaxies,
   quasars, even the centre of a supercluster. Directions on the sky (constellations,
   telescope fields, undistanced NGC entries) and Sun-centred belts cannot. */
function routable(o) {
  if (!o || !o.pos) return false;
  if (o.kind === 'constellation' || o.kind === 'deepfield' || o.kind === 'ngc') return false;
  if (o.kind === 'region' && !o.centre) return false;
  return true;
}
function routePair() {
  const a = route.from && byId[route.from], b = route.to && byId[route.to];
  return (routable(a) && routable(b) && a !== b) ? [a, b] : null;
}
function heliocentric(o) { return o.kind === 'sun' || ((o.kind === 'planet' || o.kind === 'dwarf') && o.el); }
function fmtKm(km) {
  if (km < 1e6) return Math.round(km).toLocaleString() + ' km';
  if (km < 1e9) return (km / 1e6).toFixed(1) + ' million km';
  if (km < 1e12) return (km / 1e9).toFixed(2) + ' billion km';
  if (km < 1e15) return (km / 1e12).toFixed(2) + ' trillion km';
  return (km / 1e15).toFixed(2) + ' quadrillion km';
}
function fmtDuration(s) {
  if (!isFinite(s)) return '—';
  if (s < 90) return s.toFixed(1) + ' s';
  const m = s / 60; if (m < 90) return m.toFixed(1) + ' min';
  const h = m / 60; if (h < 48) return h.toFixed(1) + ' hours';
  const dd = h / 24; if (dd < 400) return dd.toFixed(1) + ' days';
  const y = dd / 365.25;
  if (y < 1000) return y.toFixed(y < 10 ? 2 : 1) + ' years';
  if (y < 1e6) return Math.round(y).toLocaleString() + ' years';
  if (y < 1e9) return (y / 1e6).toFixed(2) + ' million years';
  if (y < 1e12) return (y / 1e9).toFixed(2) + ' billion years';
  if (y < 1e15) return (y / 1e12).toFixed(2) + ' trillion years';
  return (y / 1e15).toFixed(1) + ' quadrillion years';
}
function renderRoute() {
  const box = $('routeResult'); if (!box) return;
  const pair = routePair();
  if (!pair) { box.innerHTML = '<div class="route-hint">Choose two places that have a position in space — planets, moons, spacecraft, stars, galaxies, quasars.</div>'; return; }
  const a = pair[0], b = pair[1];
  const dAU = Math.hypot(b.pos.x - a.pos.x, b.pos.y - a.pos.y, b.pos.z - a.pos.z);
  const km = dAU * AU_KM;
  let html = '<div class="route-dist"><strong>' + esc(fmtDist(dAU)) + '</strong><span>' +
             ((dAU < 0.02 || km >= 1e15) ? '' : fmtKm(km) + ' · ') + 'right now, straight line</span></div>';
  html += '<table class="route-table">';
  VEHICLES.forEach(function (v) {
    html += '<tr><td>' + esc(v.name) + '<span>' + esc(v.note) + '</span></td><td>' + esc(fmtDuration(km / v.kmps)) + '</td></tr>';
  });
  html += '</table>';
  html += bestShot(a, b);
  const ly = dAU / LY_AU;
  if (ly > 1e8) {
    const eh = eventHorizonGly();
    html += '<div class="route-note">At these distances space itself is stretching, so these are times through a frozen universe. ' +
            (ly / 1e9 > eh
              ? '<strong>In reality nothing could ever arrive:</strong> it lies beyond the cosmic event horizon, about ' + eh.toFixed(1) + ' billion light years away — even light sent today will never reach it.'
              : 'In reality the trip would take longer, because the gap keeps growing while you travel.') + '</div>';
  }
  html += '<div class="route-note dim">Bodies keep moving — the distance updates as the clock runs.</div>';
  box.innerHTML = html;
}

/* ---- The best shot: how a real mission would actually go ---- */
const REAL_TRIPS = {
  moon: 'Apollo 11 reached lunar orbit in about 3 days (16–19 July 1969) and landed on the fourth.',
  mars: 'Perseverance took 203 days (July 2020 – February 2021), a little faster than the pure minimum-energy path.',
  venus: 'Venus Express took 153 days (November 2005 – April 2006).',
  mercury: 'MESSENGER needed 6.6 years and six planetary flybys just to slow down enough to orbit Mercury.',
  jupiter: 'Juno took five years (2011–2016), swinging past Earth once for extra speed.',
  saturn: 'Cassini took almost seven years (1997–2004), stealing speed from Venus twice, Earth and Jupiter.',
  uranus: 'Voyager 2 reached Uranus in 8.4 years (1977–1986) by slingshotting past Jupiter and Saturn.',
  neptune: 'Voyager 2 reached Neptune in 12 years (1977–1989), the only visit so far.',
  pluto: 'New Horizons flew past Pluto 9.5 years after launch (2006–2015), using a Jupiter slingshot.',
  m31: 'No need to go: Andromeda is coming to us at about 110 km/s, and the two galaxies will merge in roughly 4.5 billion years.',
  alphacenA: 'Breakthrough Starshot hopes to send gram-sized light-sail probes there at a fifth of light speed, arriving about 20 years after launch.'
};
function helioAt(o, jd) { return keplerPos(o.el, (jd - J2000) / 36525); }
function fmtDay(jd) { return dateFromJd(jd).toLocaleDateString(undefined, { day: 'numeric', month: 'short', year: 'numeric', timeZone: 'UTC' }); }
/* The next date the two bodies line up for a Hohmann transfer: the target must lead
   the departure body by pi - n_target * t_transfer at launch. */
function launchWindow(a, b) {
  const aA = a.el.a, aB = b.el.a;
  const tDays = Math.PI * Math.sqrt(Math.pow(aA + aB, 3) / (8 * MU_SUN));
  const PA = 365.25 * Math.pow(aA, 1.5), PB = 365.25 * Math.pow(aB, 1.5);
  const phi = Math.PI - TAU * tDays / PB;
  const syn = 1 / Math.abs(1 / PA - 1 / PB);
  const f = function (jd) {
    const p = helioAt(a, jd), q = helioAt(b, jd);
    let d = Math.atan2(q.y, q.x) - Math.atan2(p.y, p.x) - phi;
    return d - TAU * Math.round(d / TAU);
  };
  const step = clamp(syn / 400, 0.5, 20);
  let jd0 = simJd - step, v0 = f(jd0);
  for (let jd = simJd; jd <= simJd + syn * 1.1; jd += step) {
    const v = f(jd);
    if ((v0 <= 0) !== (v <= 0) && Math.abs(v - v0) < Math.PI) {
      let lo = jd - step, hi = jd, flo = v0;
      for (let k = 0; k < 40; k++) { const mid = (lo + hi) / 2, fm = f(mid); if ((fm <= 0) === (flo <= 0)) { lo = mid; flo = fm; } else hi = mid; }
      const launch = (lo + hi) / 2;
      return { launch: launch, arrive: launch + tDays, tDays: tDays, syn: syn };
    }
    jd0 = jd; v0 = v;
  }
  return { launch: null, tDays: tDays, syn: syn };
}
function bestShot(a, b) {
  const other = a.id === 'earth' ? b : b.id === 'earth' ? a : null;
  const real = other && REAL_TRIPS[other.id];
  let html = '';
  const orbitsSun = function (o) { return (o.kind === 'planet' || o.kind === 'dwarf') && o.el; };
  if (orbitsSun(a) && orbitsSun(b)) {
    const w = launchWindow(a, b);
    const months = w.tDays / 30.44;
    const dur = months < 24 ? months.toFixed(1) + ' months' : (w.tDays / 365.25).toFixed(1) + ' years';
    html += '<div class="route-best"><div class="route-best-title">The best shot</div>' +
      '<p>Rockets don\'t fly straight. The cheapest route is a minimum-energy transfer orbit that leaves ' + esc(a.name) +
      ' and meets ' + esc(b.name) + ' on the far side of the Sun: <strong>' + dur + '</strong> of coasting.</p>';
    if (w.launch != null) {
      const inDays = w.launch - simJd;
      const when = inDays < 45 ? 'in ' + Math.max(0, Math.round(inDays)) + ' days' : 'in ' + (inDays / 30.44).toFixed(0) + ' months';
      html += '<div class="route-window"><span>Next launch window</span><strong>' + esc(fmtDay(w.launch)) + '</strong><em>' + esc(when) + '</em></div>' +
              '<div class="route-window"><span>Arrival</span><strong>' + esc(fmtDay(w.arrive)) + '</strong><em></em></div>';
    }
    const syn = w.syn / 30.44;
    html += '<p class="dim">The planets line up like this every ' + (syn < 36 ? syn.toFixed(0) + ' months' : (w.syn / 365.25).toFixed(1) + ' years') +
            '. Launch at any other time and the trip costs far more fuel.</p>';
    if (real) html += '<p class="dim">' + esc(real) + '</p>';
    html += '</div>';
  } else if (real) {
    html += '<div class="route-best"><div class="route-best-title">How it was really done</div><p>' + esc(real) + '</p></div>';
  }
  return html;
}
let _eventHorizon = null;
/* comoving distance to the cosmic event horizon: the farthest point light sent today can ever reach */
function eventHorizonGly() {
  if (_eventHorizon == null) {
    const Ea = function (x) { const s = Math.exp(x); return Math.sqrt(COSMO.Or / (s * s * s * s) + COSMO.Om / (s * s * s) + COSMO.OL); };
    _eventHorizon = COSMO.dH * simpson(function (x) { return Math.exp(-x) / Ea(x); }, 0, 40, 8000);
  }
  return _eventHorizon;
}
function setRouteEnd(which, o) {
  route[which] = o.id;
  const inp = $(which === 'from' ? 'routeFrom' : 'routeTo'); if (inp) inp.value = o.name;
  renderRoute();
}
function openRoute(fromId, toId) {
  if (fromId) route.from = fromId;
  if (toId) route.to = toId;
  closeCompare();
  if (!route.open) route.panelWasOpen = !panel.classList.contains('hidden');
  route.open = true;
  panel.classList.add('hidden');
  $('routePanel').classList.remove('hidden');
  $('btnRoutePanel').classList.add('on');
  const f = byId[route.from], t = byId[route.to];
  $('routeFrom').value = f ? f.name : ''; $('routeTo').value = t ? t.name : '';
  renderRoute();
}
function closeRoute() {
  if (!route.open) return;
  route.open = false; $('routePanel').classList.add('hidden'); $('btnRoutePanel').classList.remove('on');
  if (route.panelWasOpen) panel.classList.remove('hidden');
}

/* ---- Size comparison: two bodies side by side, to scale ---- */
const compare = { a: 'earth', b: 'jupiter', open: false };
function sizable(o) { return o && o.radiusAU > 0 && o.kind !== 'deepfield' && o.kind !== 'constellation' && o.kind !== 'ngc' && (o.kind !== 'region' || o.centre) && o.kind !== 'distant'; }
function fmtSize(au) {
  const km = au * AU_KM;
  if (km < 1e6) return Math.round(km).toLocaleString() + ' km';
  if (km < 1e8) return (km / 1e6).toFixed(2) + ' million km';
  return fmtDist(au);
}
function sizeLabel(o) { return o.kind === 'blackhole' ? 'Shadow' : o.kind === 'region' ? 'Extent' : 'Diameter'; }
function renderCompare() {
  const box = $('compareResult'); if (!box) return;
  const A = byId[compare.a], B = byId[compare.b];
  if (A && A === B) { box.innerHTML = '<div class="route-hint">Pick two different things to compare.</div>'; return; }
  if (!sizable(A) || !sizable(B)) { box.innerHTML = '<div class="route-hint">Pick two things with a size — planets, moons, stars, black holes, galaxies, even the Local Group.</div>'; return; }
  const big = A.radiusAU >= B.radiusAU ? A : B, small = big === A ? B : A;
  const ratio = big.radiusAU / small.radiusAU, vol = Math.pow(ratio, 3);
  const fmtN = function (v) { return v >= 1e12 ? (v / 1e12).toFixed(1) + ' trillion' : v >= 1e9 ? (v / 1e9).toFixed(1) + ' billion' : v >= 1e6 ? (v / 1e6).toFixed(1) + ' million' : v >= 100 ? Math.round(v).toLocaleString() : v.toFixed(v >= 10 ? 1 : 2); };
  const cw = Math.max(240, box.clientWidth || 306), ch = Math.round(cw * 0.56), dpr = Math.min(window.devicePixelRatio || 1, 2);
  let html = '<canvas id="compareCanvas" width="' + Math.round(cw * dpr) + '" height="' + Math.round(ch * dpr) + '" style="height:' + ch + 'px"></canvas>';
  html += '<table class="route-table"><tr><td>' + esc(A.name) + '<span>' + esc(sizeLabel(A)) + '</span></td><td>' + esc(fmtSize(A.radiusAU * 2)) + '</td></tr>' +
          '<tr><td>' + esc(B.name) + '<span>' + esc(sizeLabel(B)) + '</span></td><td>' + esc(fmtSize(B.radiusAU * 2)) + '</td></tr></table>';
  html += '<div class="route-best"><div class="route-best-title">Side by side</div>' +
          (ratio < 1.05
            ? '<p><strong>' + esc(big.name) + '</strong> and <strong>' + esc(small.name) + '</strong> are about the same size.</p>'
            : '<p><strong>' + esc(big.name) + '</strong> is <strong>' + fmtN(ratio) + '×</strong> wider than ' + esc(small.name) + '.' +
              (vol < 2 ? ' It has about <strong>' + vol.toFixed(1) + '×</strong> the volume.' : ratio < 1e5 ? ' About <strong>' + fmtN(vol) + '</strong> ' + esc(/s$/i.test(small.name) ? small.name : small.name + 's') + ' would fit inside it.' : '') + '</p>');
  if (byId.earth && A !== byId.earth && B !== byId.earth) {
    const eA = A.radiusAU / byId.earth.radiusAU, eB = B.radiusAU / byId.earth.radiusAU;
    html += '<p class="dim">For scale: ' + esc(A.name) + ' is ' + fmtN(eA) + '× Earth, ' + esc(B.name) + ' is ' + fmtN(eB) + '× Earth.</p>';
  }
  if (ratio > 40) html += '<p class="dim">At this scale ' + esc(small.name) + ' is a dot, so it is drawn with a magnified inset.</p>';
  html += '</div>';
  box.innerHTML = html;
  drawCompare($('compareCanvas'), A, B, big, small, ratio, cw, ch, dpr);
}
function bodyLook(o) { return { color: o.color, color2: o.color2 || o.color, rings: o.rings ? o.rings.outerKm / (o.radiusAU * AU_KM) : 0, bh: o.kind === 'blackhole', region: o.kind === 'region', star: o.kind === 'star' || o.kind === 'sun' }; }
function drawBall(g2, x, y, r, look) {
  if (look.region) {
    const gr = g2.createRadialGradient(x, y, r * 0.5, x, y, r); gr.addColorStop(0, hexA(look.color, 0.06)); gr.addColorStop(0.92, hexA(look.color, 0.35)); gr.addColorStop(1, hexA(look.color, 0));
    g2.fillStyle = gr; g2.beginPath(); g2.arc(x, y, r, 0, TAU); g2.fill(); return;
  }
  if (look.bh) {
    const gr = g2.createRadialGradient(x, y, r, x, y, r * 1.35); gr.addColorStop(0, 'rgba(255,190,110,0.9)'); gr.addColorStop(1, 'rgba(255,190,110,0)');
    g2.fillStyle = gr; g2.beginPath(); g2.arc(x, y, r * 1.35, 0, TAU); g2.fill();
    g2.fillStyle = '#000'; g2.beginPath(); g2.arc(x, y, r, 0, TAU); g2.fill(); return;
  }
  if (look.rings > 1) {
    g2.save(); g2.translate(x, y); g2.scale(1, 0.28);
    g2.strokeStyle = hexA(look.color2, 0.55); g2.lineWidth = r * (look.rings - 1.28) ; g2.beginPath(); g2.arc(0, 0, r * (1.28 + look.rings) / 2, 0, TAU); g2.stroke(); g2.restore();
  }
  const gr = g2.createRadialGradient(x - r * 0.35, y - r * 0.35, r * 0.05, x, y, r);
  gr.addColorStop(0, look.star ? '#ffffff' : look.color); gr.addColorStop(0.6, look.color); gr.addColorStop(1, look.star ? look.color2 : 'rgba(0,0,0,0.85)');
  g2.fillStyle = gr; g2.beginPath(); g2.arc(x, y, r, 0, TAU); g2.fill();
  if (look.star) { const gl = g2.createRadialGradient(x, y, r, x, y, r * 1.5); gl.addColorStop(0, hexA(look.color, 0.5)); gl.addColorStop(1, hexA(look.color, 0)); g2.fillStyle = gl; g2.beginPath(); g2.arc(x, y, r * 1.5, 0, TAU); g2.fill(); }
}
function drawCompare(cv, A, B, big, small, ratio, W2, H2, dpr) {
  if (!cv) return;
  const g2 = cv.getContext('2d');
  g2.setTransform(dpr, 0, 0, dpr, 0, 0);
  g2.clearRect(0, 0, W2, H2);
  g2.fillStyle = '#07070b'; g2.fillRect(0, 0, W2, H2);
  const lookA = bodyLook(A), lookB = bodyLook(B);
  const bigLook = A === big ? lookA : lookB, smallLook = A === big ? lookB : lookA;
  let rBig = Math.min(H2 * 0.36, W2 * 0.26);
  if (bigLook.rings > 1) rBig = Math.min(rBig, W2 * 0.26 / bigLook.rings);
  const rSmall = Math.max(rBig / ratio, 0.6);
  const cx1 = W2 * 0.33, cx2 = W2 * 0.72, cy = H2 * 0.5;
  const leftIsA = A === big;
  drawBall(g2, cx1, cy, rBig, leftIsA ? lookA : lookB);
  drawBall(g2, cx2, cy, rSmall, leftIsA ? lookB : lookA);
  if (ratio > 40) {                                          /* inset for the speck */
    const ir = Math.round(H2 * 0.13) / Math.max(1, smallLook.rings), ix = cx2, iy = cy - H2 * 0.27;
    g2.strokeStyle = 'rgba(255,255,255,0.25)'; g2.lineWidth = 1; g2.beginPath(); g2.arc(ix, iy, ir + 6, 0, TAU); g2.stroke();
    g2.beginPath(); g2.moveTo(ix, iy + ir + 6); g2.lineTo(cx2, cy - rSmall - 2); g2.stroke();
    drawBall(g2, ix, iy, ir, leftIsA ? lookB : lookA);
    g2.fillStyle = 'rgba(255,255,255,0.55)'; g2.font = '500 12px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'; g2.textAlign = 'center';
    g2.fillText('magnified ' + Math.round(ir * ratio / rBig).toLocaleString() + '×', ix, iy + ir + 22);
  }
  g2.fillStyle = '#e8e8f0'; g2.font = '600 12.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'; g2.textAlign = 'center';
  g2.fillText(big.name, cx1, H2 - 12); g2.fillText(small.name, cx2, H2 - 12);
}
function openCompare(aId, bId) {
  if (aId) compare.a = aId;
  if (bId) compare.b = bId;
  closeRoute();
  if (!compare.open) compare.panelWasOpen = !panel.classList.contains('hidden');
  compare.open = true;
  panel.classList.add('hidden');
  $('comparePanel').classList.remove('hidden');
  $('btnComparePanel').classList.add('on');
  const A = byId[compare.a], B = byId[compare.b];
  $('compareA').value = A ? A.name : ''; $('compareB').value = B ? B.name : '';
  renderCompare();
}
function closeCompare() {
  if (!compare.open) return;
  compare.open = false; $('comparePanel').classList.add('hidden'); $('btnComparePanel').classList.remove('on');
  if (compare.panelWasOpen) panel.classList.remove('hidden');
}
/* a small search picker bound to one input */
function attachPicker(inputId, listId, which, accept, onPick) {
  const inp = $(inputId), list = $(listId);
  accept = accept || routable;
  function run() {
    const q = inp.value.trim().toLowerCase();
    if (!q) { list.classList.remove('show'); return; }
    let rows = objects.filter(function (o) { return accept(o) && o.name.toLowerCase().indexOf(q) >= 0; }).slice(0, 8);
    if (rows.length < 8) searchBulk(q, rows, 8);
    rows = rows.filter(function (r) { return r.bulkRow ? !(accept === sizable && /^(sat|ngc)/.test(r.id)) : accept(r); });
    list.innerHTML = '';
    rows.forEach(function (r) {
      const div = document.createElement('div'); div.className = 'search-row';
      div.innerHTML = '<span class="dot" style="background:' + esc(r.color) + '"></span><span class="nm">' + esc(r.name) + '</span><span class="ty">' + esc(r.type) + '</span>';
      div.addEventListener('click', function () { const o = r.bulkRow ? r.promote() : r; if (!accept(o)) { toast('No size known for ' + o.name); return; } (onPick || setRouteEnd)(which, o); list.classList.remove('show'); });
      list.appendChild(div);
    });
    list.classList.toggle('show', rows.length > 0);
  }
  inp.addEventListener('input', run);
  inp.addEventListener('focus', function () { inp.select(); run(); });
  inp.addEventListener('keydown', function (e) { if (e.key === 'Escape') { list.classList.remove('show'); inp.blur(); } });
}
/* the route drawn on the map */
function drawRoute() {
  if (!route.open) return;
  const pair = routePair(); if (!pair) return;
  const pa = project(pair[0].pos), pb = project(pair[1].pos);
  if (!pa || !pb) return;
  ctx.save();
  ctx.strokeStyle = 'rgba(124,140,255,0.85)'; ctx.lineWidth = 1.5; ctx.setLineDash([6, 5]);
  ctx.beginPath(); ctx.moveTo(pa.x, pa.y); ctx.lineTo(pb.x, pb.y); ctx.stroke();
  ctx.setLineDash([]);
  [pa, pb].forEach(function (p) { ctx.fillStyle = '#7c8cff'; ctx.beginPath(); ctx.arc(p.x, p.y, 4, 0, TAU); ctx.fill(); ctx.strokeStyle = '#fff'; ctx.lineWidth = 1; ctx.stroke(); });
  const mx = (pa.x + pb.x) / 2, my = (pa.y + pb.y) / 2;
  const dAU = Math.hypot(pair[1].pos.x - pair[0].pos.x, pair[1].pos.y - pair[0].pos.y, pair[1].pos.z - pair[0].pos.z);
  const txt = fmtDist(dAU);
  ctx.font = '600 11px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
  const w = ctx.measureText(txt).width + 14;
  ctx.fillStyle = 'rgba(16,16,22,0.9)'; ctx.strokeStyle = 'rgba(124,140,255,0.6)';
  ctx.beginPath(); ctx.roundRect ? ctx.roundRect(mx - w / 2, my - 10, w, 20, 10) : ctx.rect(mx - w / 2, my - 10, w, 20); ctx.fill(); ctx.stroke();
  ctx.fillStyle = '#c3caff'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(txt, mx, my);
  ctx.restore();
}

/* ---- search ---- */
const searchInput = $('searchInput'), searchResults = $('searchResults');
let searchSel = -1, searchList = [];

function runSearch() {
  const q = searchInput.value.trim().toLowerCase();
  if (!q) { searchResults.classList.remove('show'); searchList = []; searchSel = -1; return; }
  /* exact names and ids first, then names that start with the query, then anything
     containing it — so "leo" is Leo before Galileo and "m31" is Andromeda */
  const qk = q.replace(/\s+/g, '');
  const scored = [];
  for (let i = 0; i < objects.length; i++) {
    const o = objects[i], nm = o.name.toLowerCase();
    let s = -1;
    if (nm === q || o.id.toLowerCase() === qk) s = 0;
    else if (nm.indexOf(q) === 0) s = 1;
    else if (nm.indexOf(q) >= 0) s = 2;
    else if ((o.type || '').toLowerCase().indexOf(q) >= 0) s = 3;
    const al = o.data && o.data.aliases;
    if (al) for (let k = 0; k < al.length; k++) {
      const v = al[k].toLowerCase();
      if (v === q || v.replace(/\s+/g, '') === qk) { s = 0; break; }
      if (v.indexOf(q) === 0 && (s < 0 || s > 1)) s = 1;
    }
    if (s >= 0) scored.push({ o: o, s: s, p: o.prio || 0 });
  }
  scored.sort(function (x, y) { return x.s - y.s || y.p - x.p; });
  searchList = scored.slice(0, 12).map(function (x) { return x.o; });
  if (searchList.length < 12) searchBulk(q, searchList, 12);
  searchResults.innerHTML = '';
  if (!searchList.length) {
    searchResults.innerHTML = '<div class="search-empty">Nothing found for "' + esc(q) + '"</div>';
  } else {
    searchList.forEach(function (o, i) {
      const row = document.createElement('div');
      row.className = 'search-row' + (i === searchSel ? ' sel' : '');
      row.innerHTML = '<span class="dot" style="background:' + esc(o.color) + '"></span>' +
                      '<span class="nm">' + esc(o.name) + '</span><span class="ty">' + esc(o.type) + '</span>';
      row.addEventListener('click', function () { goSearch(o); });
      searchResults.appendChild(row);
    });
  }
  searchResults.classList.add('show');
}
function goSearch(o) {
  if (o.bulkRow) o = o.promote();   /* only now does it enter the scene */
  select(o.id); flyTo(o.id);
  searchResults.classList.remove('show');
  searchInput.blur();
}
searchInput.addEventListener('input', function () { searchSel = -1; runSearch(); });
searchInput.addEventListener('focus', runSearch);
searchInput.addEventListener('keydown', function (e) {
  if (e.key === 'ArrowDown') { searchSel = Math.min(searchSel + 1, searchList.length - 1); runSearch(); e.preventDefault(); }
  else if (e.key === 'ArrowUp') { searchSel = Math.max(searchSel - 1, 0); runSearch(); e.preventDefault(); }
  else if (e.key === 'Enter') { const o = searchList[Math.max(searchSel, 0)]; if (o) goSearch(o); }
  else if (e.key === 'Escape') { searchResults.classList.remove('show'); searchInput.blur(); }
});
document.addEventListener('click', function (e) {
  if (!e.target.closest('.search-wrap')) searchResults.classList.remove('show');
  if (!e.target.closest('#settingsPopover') && e.target.closest('#btnSettings') === null) $('settingsPopover').classList.add('hidden');
  if (!e.target.closest('#eventsPopover') && e.target.closest('#btnOnThisDay') === null) $('eventsPopover').classList.add('hidden');
  if (!e.target.closest('#aboutPopover') && e.target.closest('#btnMore') === null) $('aboutPopover').classList.add('hidden');
});

/* ---- trip planner wiring ---- */
$('btnRoutePanel').addEventListener('click', function () { if (route.open) closeRoute(); else { closeCompare(); openRoute(null, route.to || selectedId); } });
$('routeClose').addEventListener('click', closeRoute);
$('routeSwap').addEventListener('click', function () { const t = route.from; route.from = route.to; route.to = t; openRoute(); });
attachPicker('routeFrom', 'routeFromList', 'from');
attachPicker('routeTo', 'routeToList', 'to');
/* size comparison wiring */
function setCompareEnd(which, o) { compare[which] = o.id; $(which === 'a' ? 'compareA' : 'compareB').value = o.name; renderCompare(); }
attachPicker('compareA', 'compareAList', 'a', sizable, setCompareEnd);
attachPicker('compareB', 'compareBList', 'b', sizable, setCompareEnd);
$('btnComparePanel').addEventListener('click', function () {
  if (compare.open) closeCompare();
  else if (selectedId && sizable(byId[selectedId])) openCompare(selectedId, selectedId === compare.b ? (selectedId === 'earth' ? 'jupiter' : 'earth') : null);
  else openCompare();
});
$('compareClose').addEventListener('click', closeCompare);
$('compareSwap').addEventListener('click', function () { const t = compare.a; compare.a = compare.b; compare.b = t; openCompare(); });
document.querySelectorAll('#compareQuick button').forEach(function (b) {
  b.addEventListener('click', function () { setCompareEnd('b', byId[b.getAttribute('data-id')]); });
});
$('btnCompare').addEventListener('click', function () { if (selectedId && sizable(byId[selectedId])) openCompare(selectedId, selectedId === 'earth' ? 'jupiter' : 'earth'); });
$('btnRoute').addEventListener('click', function () { if (selectedId) openRoute(selectedId === 'earth' ? 'sun' : 'earth', selectedId); });
document.addEventListener('click', function (e) {
  if (!e.target.closest('.picker')) document.querySelectorAll('.picker-list').forEach(function (l) { l.classList.remove('show'); });
});
/* checkbox captions carry the real counts */
if (bulkStars) $('lblAllStars').textContent = 'Deep star field — all ' + bulkStars.n.toLocaleString() + ' stars (default ' + bulkStars.defaultN.toLocaleString() + ')';
if (bulkGal) $('lblAllGalaxies').textContent = 'All ' + bulkGal.n.toLocaleString() + ' galaxies (default ' + bulkGal.defaultN.toLocaleString() + ')';
if (bulkSbdb) $('lblAllSmall').textContent = 'All ' + bulkSbdb.n.toLocaleString() + ' asteroids & comets (default ' + bulkSbdb.nImp.toLocaleString() + ')';

/* ---- history ---- */
const history = [];
let historyIdx = -1;
let suppressHistory = false;
function pushHistory(id) {
  if (suppressHistory) return;   /* navigating history is not itself history */
  if (history[historyIdx] === id) return;
  history.splice(historyIdx + 1);
  history.push(id);
  historyIdx = history.length - 1;
  updateHistoryButtons();
}
function updateHistoryButtons() {
  $('btnBack').disabled = historyIdx <= 0;
  $('btnForward').disabled = historyIdx >= history.length - 1;
}
$('btnBack').addEventListener('click', function () {
  if (historyIdx <= 0) return;
  historyIdx--; const id = history[historyIdx];
  select(id); flyToNoHistory(id); updateHistoryButtons();
});
$('btnForward').addEventListener('click', function () {
  if (historyIdx >= history.length - 1) return;
  historyIdx++; const id = history[historyIdx];
  select(id); flyToNoHistory(id); updateHistoryButtons();
});
function flyToNoHistory(id) {
  suppressHistory = true;
  try { flyTo(id); } finally { suppressHistory = false; }
  updateHistoryButtons();
}

/* ---- home ---- */
function goHome() {
  cam.follow = null; cam.flight = null; selectedId = null; focusCategory = null;
  $('infoPanel').classList.add('hidden');
  panelDetail.classList.add('hidden');
  panelRoot.classList.remove('hidden');
  cam.gx = cam.gy = cam.gz = 0;
  cam.distGoal = 48; cam.pitchGoal = 0.42; cam.fovGoal = FOV_DEFAULT;
  cam.rollGoal = 0; cam.pendYaw = cam.pendPitch = 0;
  updateFocusChip();
}
$('btnHome').addEventListener('click', goHome);
$('brandHome').addEventListener('click', goHome);

/* ---- zoom ---- */
$('zoomIn').addEventListener('click', function () { zoomBy(1 / 1.7); });
$('zoomOut').addEventListener('click', function () { zoomBy(1.7); });

/* ---- settings ---- */
const optMap = { optOrbits: 'orbits', optLabels: 'labels', optMoons: 'moons', optBelts: 'belts',
                 optOort: 'oort', optStars: 'stars', optMilkyWay: 'milkyway', optProbes: 'probes',
                 optRealSize: 'realSize', optConstellations: 'constellations', optDeepFields: 'deepfields', optNgc: 'ngc',
                 optMainBelt: 'mainbelt', optNeo: 'neo', optTrojans: 'trojans', optTno: 'tno', optComets: 'comets',
                 optGalaxies: 'galaxies', optSatellites: 'satellites',
                 optAllStars: 'allStars', optAllGalaxies: 'allGalaxies', optAllSmall: 'allSmallBodies', optSurfaces: 'surfaces', optFreeRotate: 'freeRotate' };
Object.keys(optMap).forEach(function (elId) {
  const el = $(elId);
  el.checked = opts[optMap[elId]];
  el.addEventListener('change', function () { opts[optMap[elId]] = el.checked; });
});
document.querySelectorAll('#segStructures button').forEach(function (b) {
  b.classList.toggle('on', b.getAttribute('data-v') === opts.structures);
  b.addEventListener('click', function () {
    opts.structures = b.getAttribute('data-v');
    document.querySelectorAll('#segStructures button').forEach(function (x) { x.classList.toggle('on', x === b); });
  });
});
$('btnSettings').addEventListener('click', function (e) {
  e.stopPropagation();
  const p = $('settingsPopover');
  p.classList.toggle('hidden');
  $('eventsPopover').classList.add('hidden');
  $('aboutPopover').classList.add('hidden');
});

/* ---- on this day ---- */
$('btnOnThisDay').addEventListener('click', function (e) {
  e.stopPropagation();
  const p = $('eventsPopover');
  const wasHidden = p.classList.contains('hidden');
  $('settingsPopover').classList.add('hidden');
  $('aboutPopover').classList.add('hidden');
  p.classList.toggle('hidden');
  if (!wasHidden) return;
  const d = dateFromJd(simJd);
  const key = String(d.getUTCMonth() + 1).padStart(2, '0') + '-' + String(d.getUTCDate()).padStart(2, '0');
  $('eventsTitle').textContent = 'On ' + d.toLocaleDateString(undefined, { day: 'numeric', month: 'long', timeZone: 'UTC' });
  const ul = $('eventList');
  ul.innerHTML = '';
  const list = EVENTS[key];
  if (!list) {
    ul.innerHTML = '<li style="color:var(--muted)">No milestone recorded for this date — try 20 July, 5 September or 14 July.</li>';
  } else {
    list.forEach(function (ev) {
      const li = document.createElement('li');
      li.innerHTML = '<span class="event-year">' + ev.y + '</span>' + ev.t;
      ul.appendChild(li);
    });
  }
});
$('btnMore').addEventListener('click', function (e) {
  e.stopPropagation();
  $('settingsPopover').classList.add('hidden');
  $('eventsPopover').classList.add('hidden');
  $('aboutPopover').classList.toggle('hidden');
});

/* ---- time ---- */
const SPEEDS = [
  { label: '1×', v: 1 }, { label: '1 h/s', v: 3600 }, { label: '1 d/s', v: 86400 },
  { label: '1 wk/s', v: 604800 }, { label: '1 mo/s', v: 2629800 }, { label: '1 yr/s', v: 31557600 }
];
let speedIdx = 2, playing = true;

SPEEDS.forEach(function (s, i) {
  const b = document.createElement('button');
  b.className = 'speed-btn' + (i === speedIdx ? ' on' : '');
  b.textContent = s.label;
  b.addEventListener('click', function () {
    speedIdx = i;
    playing = true;
    document.querySelectorAll('.speed-btn').forEach(function (x, j) { x.classList.toggle('on', j === i); });
    updatePlayIcon();
  });
  $('timeSpeeds').appendChild(b);
});
function updatePlayIcon() {
  $('iconPlay').classList.toggle('hidden', playing);
  $('iconPause').classList.toggle('hidden', !playing);
}
$('btnPlay').addEventListener('click', function () { playing = !playing; updatePlayIcon(); });
$('btnNow').addEventListener('click', function () { simJd = jdFromDate(new Date()); });
$('btnTime').addEventListener('click', function () {
  $('timeBar').classList.toggle('hidden');
  $('btnTime').classList.toggle('on', !$('timeBar').classList.contains('hidden'));
});
$('btnTime').classList.add('on');

/* ---- toast + hint ---- */
let toastTimer = null;
function toast(msg) {
  const t = $('toast');
  t.textContent = msg;
  t.classList.add('show');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(function () { t.classList.remove('show'); }, 2200);
}
setTimeout(function () { $('hint').classList.add('show'); }, 900);
setTimeout(function () { $('hint').classList.remove('show'); }, 6500);

/* ============================================================
   9. Main loop
   ============================================================ */
let last = performance.now();
let statTick = 0;

function frame(now) {
  const dt = Math.min((now - last) / 1000, 0.1);
  last = now;

  if (deep.on) deepTick(dt);
  else if (playing) simJd += (SPEEDS[speedIdx].v * dt) / 86400;

  updatePositions(simJd);
  if (deep.on) deepWorld();

  /* a slow cinematic drift until the viewer takes over */
  idle += dt;
  if (idle > 6 && !dragging && !cam.flight && cam.follow === null && !deep.on) {
    const drift = dt * 0.012 * (cam.fov / FOV_DEFAULT);
    if (opts.freeRotate && cam.fovGoal >= FOV_DEFAULT * 0.999) cam.pendYaw += drift; else cam.yawGoal += drift;
  }

  if (deep.on && !deep.takeover) deepCamera(dt);
  updateCamera(dt);
  render(dt);
  if (deep.on) { deepOverlay(); deepHud(); }

  if (route.open && statTick % 30 === 0) renderRoute();
  /* distances drift as the clock runs, so refresh them a few times a second */
  if (selectedId && byId[selectedId] && ++statTick % 20 === 0) renderStats(byId[selectedId]);
  $('timeDate').textContent = fmtDate(simJd);
  $('scaleReadout').textContent = cam.fov < FOV_DEFAULT * 0.98 ? 'View: ' + fmtAngle(cam.fov) + ' across' : 'View: ' + fmtDist(cam.dist);

  requestAnimationFrame(frame);
}

/* promote the famous small bodies so the Asteroids and Comets categories are populated */
if (bulkSbdb && typeof FEATURED_ASTEROIDS !== 'undefined') {
  const wantA = {}, wantC = {};
  FEATURED_ASTEROIDS.forEach(function (nm) { wantA[nm.toLowerCase()] = true; });
  FEATURED_COMETS.forEach(function (nm) { wantC[nm.toLowerCase()] = true; });
  const takenA = {}, takenC = {};
  for (let k = 0; k < bulkSbdb.n; k++) {
    const comet = bulkSbdb.cls[k] >= 4;
    const want = comet ? wantC : wantA, taken = comet ? takenC : takenA;
    const nm = want[bulkSbdb.lkey[k]] ? bulkSbdb.lkey[k] : (want[bulkSbdb.lname[k]] ? bulkSbdb.lname[k] : null);
    if (nm && !taken[nm]) { bulkSbdbObject(k); taken[nm] = true; }
  }
}
if (bulkSat && typeof FEATURED_SATELLITES !== 'undefined') {
  const taken = {}; objects.forEach(function (o) { taken[o.name.toLowerCase()] = true; });
  FEATURED_SATELLITES.forEach(function (pref) {
    const p = pref.toLowerCase();
    for (let k = 0; k < bulkSat.n; k++) {
      if (bulkSat.lname[k].indexOf(p) === 0 && !taken[bulkSat.lname[k]]) { bulkSatObject(k); taken[bulkSat.lname[k]] = true; break; }
    }
  });
}
resize();
updatePositions(simJd);
updateCamera(0.016);
updatePlayIcon();
updateHistoryButtons();
requestAnimationFrame(frame);

/* expose a little of the internals for tinkering from the console */

/* ============================================================
   Deep time — the future of the universe on a log clock
   ============================================================ */
const DEEP_SHOTS_FN = {};
const deepStore = {};
function fmtYears(Y) {
  if (Y < 1e4) return Math.round(Y).toLocaleString() + (Math.round(Y) === 1 ? ' year' : ' years');
  if (Y < 1e6) return Math.round(Y / 1e3).toLocaleString() + ' thousand years';
  if (Y < 1e9) return (Y / 1e6).toFixed(Y < 1e7 ? 1 : 0) + ' million years';
  if (Y < 1e12) return (Y / 1e9).toFixed(Y < 1e10 ? 2 : 0) + ' billion years';
  if (Y < 1e15) return (Y / 1e12).toFixed(Y < 1e13 ? 1 : 0) + ' trillion years';
  return null;                                        /* shown as a power of ten */
}
/* how long each catalogue star has left, from its colour and brightness */
function massFromBV(bv) {
  const T = [[-0.33, 20], [-0.2, 7], [-0.1, 4], [0, 2.5], [0.15, 1.9], [0.3, 1.5], [0.45, 1.25], [0.6, 1.05], [0.65, 1.0], [0.8, 0.9], [1.0, 0.78], [1.2, 0.68], [1.4, 0.6], [1.6, 0.45], [1.8, 0.3], [2.0, 0.18]];
  if (bv <= T[0][0]) return T[0][1];
  for (let i = 1; i < T.length; i++) if (bv <= T[i][0]) { const f = (bv - T[i - 1][0]) / (T[i][0] - T[i - 1][0]); return T[i - 1][1] + f * (T[i][1] - T[i - 1][1]); }
  return 0.15;
}
function msLife(m) { return m >= 0.5 ? 1e10 * Math.pow(m, -2.5) : 5.7e10 * Math.pow(0.5 / m, 3.2); }
function starLives() {
  const S = bulkStars;
  if (S.life) return S.life;
  const L = new Float32Array(S.n);
  for (let i = 0; i < S.n; i++) {
    const dpc = Math.max(Math.hypot(S.x[i], S.y[i], S.z[i]) / PC_AU, 1e-3);
    const M = S.mag[i] * 0.01 - 5 * Math.log10(dpc) + 5, bv = S.ci[i] * 0.001;
    let life;
    if (M < -4) life = 1e6;                                  /* supergiants */
    else if (bv > 0.75 && M < 1.2) life = 3e8;               /* giants: already old */
    else life = msLife(massFromBV(bv)) * 0.5;                /* on average half way through */
    L[i] = life;
  }
  S.life = L;
  return L;
}
function curatedLife(o) {
  if (o.id === 'betelgeuse') return 1e5;
  if (o.category === 'remnants') return 1e16;
  const d = o.data || {}, r = d.rSun || 1, sp = (d.spec || 'G').charAt(0);
  if (r > 300) return 3e5;
  if (r > 50) return 2e6;
  if (r > 10) return 3e8;
  return ({ O: 5e6, B: 1e8, A: 1e9, F: 4e9, G: 6e9, K: 2e10, M: 1e12 }[sp] || 6e9);
}
function sunAt(Y) {
  const T = SUN_TRACK;
  if (Y >= T[T.length - 1][0]) return { r: T[T.length - 1][1], col: T[T.length - 1][2] };
  for (let i = 1; i < T.length; i++) if (Y <= T[i][0]) {
    const f = (Y - T[i - 1][0]) / (T[i][0] - T[i - 1][0]);
    const r = Math.exp(lerp(Math.log(T[i - 1][1]), Math.log(T[i][1]), f));
    const c0 = parseInt(T[i - 1][2].slice(1), 16), c1 = parseInt(T[i][2].slice(1), 16);
    const ch = function (s) { return Math.round(lerp(c0 >> s & 255, c1 >> s & 255, f)); };
    return { r: r, col: '#' + [ch(16), ch(8), ch(0)].map(function (v) { return v.toString(16).padStart(2, '0'); }).join('') };
  }
  return { r: 1, col: '#ffd977' };
}
function sunMaxRUpTo(Y) {
  let m = sunAt(Y).r;
  SUN_TRACK.forEach(function (k) { if (k[0] <= Y) m = Math.max(m, k[1]); });
  return m;
}
/* planets, moons and comets that no longer exist or move too fast to see */
function deepGone(o) {
  const Y = deep.Y;
  if (Y > 1e15 && (o.kind === 'planet' || o.kind === 'dwarf' || o.kind === 'moon' || o.kind === 'exoplanet')) return true;   /* flung away by passing stars */
  const host = o.kind === 'moon' ? byId[o.parent] : o;
  if (host && host.el && (host.kind === 'planet') && host.el.a * 215 < deep.sunMaxR * 1.02) return true;                    /* swallowed by the red giant */
  return false;
}
function deepHidden(o) {
  const Y = deep.Y;
  if (o.kind === 'probe' || o.category === 'satellites') return true;
  if (o.kind === 'planet' || o.kind === 'dwarf' || o.kind === 'moon' || o.kind === 'exoplanet') {
    if (o.id === 'earth' && deep.shot === 'earth') return false;   /* the Earth chapter shows Earth itself */
    return deepGone(o) || !!o.deepBlur;
  }
  if (o.kind === 'sun') return Y > 1e16;
  if (o.kind === 'star') { if (o.deepLife == null) o.deepLife = curatedLife(o); return Y > o.deepLife; }
  if (o.kind === 'blackhole') { if (o.deepEvap == null) o.deepEvap = 2.1e67 * Math.pow((o.data && o.data.massSuns) || 10, 3); return Y > o.deepEvap; }
  if (o.kind === 'neutron') return Y > 1e16;
  if (o.kind === 'deepsky') {
    if (o.category === 'galaxies') { if (o.id === 'm31') return deep.merge > 0.9; return o.lg ? Y > 8e9 : deep.horizonFade < 0.4; }   /* the satellites end up inside Milkomeda */
    return Y > 1e8;                                         /* nebulae disperse, clusters drift apart */
  }
  if (o.kind === 'distant') return deep.horizonFade < 0.4;
  return false;
}
/* ---- the clock: ten times faster every second, pausing at each chapter ---- */
function deepRate(L) {
  const zones = [[7.0, 8.36, 0.3], [9.15, 9.65, 0.22], [9.74, 9.879, 0.07], [10.3, 11.18, 0.45],
                 [20.5, 33.6, 3], [41, 66.5, 3.5], [68.5, 86.6, 3.5], [88, 99.8, 3.5]];
  for (let i = 0; i < zones.length; i++) if (L >= zones[i][0] && L < zones[i][1]) return zones[i][2];
  return 1;
}
function deepSetLog(L) {
  deep.logY = clamp(L, 0, 100.3);
  deep.Y = Math.pow(10, deep.logY);
}
function deepTick(dt) {
  deep.clock += dt;
  let rate = 0;
  if (deep.playing) {
    if (deep.hold > 0) { deep.hold -= dt; rate = 0.012; }
    else rate = deepRate(deep.logY);
    const next = DEEP_TIME[deep.chapter + 1];
    let L = deep.logY + rate * dt;
    if (deep.hold > 0 && next) L = Math.min(L, next.log - 1e-4);   /* the slow drift during a pause never reaches the next chapter */
    if (next && L >= next.log) {                             /* arrive at a chapter and dwell on it */
      L = next.log; deep.chapter++;
      const words = next.text.split(' ').length + next.title.split(' ').length;
      deep.hold = clamp(1.6 + words * 0.13, 3, 6.5);
      deepShowChapter();
    }
    deepSetLog(L);
    if (deep.chapter >= DEEP_TIME.length - 1 && deep.hold <= 0) { deep.playing = false; deepEnd(true); }
  }
  deep.yrsPerSec = deep.Y * Math.LN10 * rate;
  /* the planets really move for the first thousand years; after that they are rings */
  simJd = deep.jd0 + Math.min(deep.Y, 1000) * 365.25;
}
function deepWorld() {
  const Y = deep.Y;
  const yrsPerFrame = deep.yrsPerSec / 60;
  objects.forEach(function (o) {
    if (o.kind !== 'planet' && o.kind !== 'dwarf' && o.kind !== 'moon') return;
    if (o.deepBlur) return;
    const P = o.el ? Math.pow(o.el.a, 1.5) : (o.periodDays ? o.periodDays / 365.25 : 1e9);
    if (yrsPerFrame > P / 6 || Y > 1000) o.deepBlur = true;   /* sticky: once a ring, stays a ring */
  });
  const s = sunAt(Y), sun = byId.sun;
  sun.radiusAU = s.r * R_SUN_AU; sun.color = s.col;
  deep.sunMaxR = sunMaxRUpTo(Y);
  const e = byId.earth;
  if (Y > 1e9) { e.dry = true; e.color = '#a8805c'; e.color2 = '#5a3c24'; }
  else { e.dry = false; e.color = deepStore.earthColor; e.color2 = deepStore.earthColor2; }
  deep.horizonFade = Y < 2e10 ? 1 : clamp(1 - (Math.log10(Y) - 10.3) / 0.88, 0, 1);
  deep.galaxyFade = Y < 1e12 ? 1 : Y < 1e14 ? lerp(1, 0.15, (Math.log10(Y) - 12) / 2) : clamp(0.15 * (1 - (Math.log10(Y) - 14) / 6), 0, 0.15);
  deep.merge = clamp((Y - 4.9e9) / 1.6e9, 0, 1);
  /* Andromeda falls toward us */
  const m31 = byId.m31;
  if (m31 && MW) {
    const s2 = clamp(1 - Math.pow(Y / 4.9e9, 2), 0, 1), c = MW.centre, p0 = deepStore.m31pos;   /* first close pass about 4.5 billion years out, merged by about 6.5 */
    m31.pos = { x: c.x + (p0.x - c.x) * s2, y: c.y + (p0.y - c.y) * s2, z: c.z + (p0.z - c.z) * s2 };
  }
}
/* ---- the guided camera ---- */
function deepShot(name) {
  const Y = deep.Y, t = deep.clock;
  const ang = function (yaw, pitch) { return { yaw: yaw, pitch: pitch, roll: 0 }; };
  if (name === 'solar') return Object.assign({ pos: { x: 0, y: 0, z: 0 }, dist: 48 }, ang(-1.05 + 0.04 * t, 0.5));
  if (name === 'solarWide') return Object.assign({ pos: { x: 0, y: 0, z: 0 }, dist: 7.5 }, ang(-1.05 + 0.05 * t, 0.55));
  if (name === 'earth') { const e = byId.earth; const v = sunlitAngles(e) || { yaw: 0, pitch: 0.3 }; return Object.assign({ pos: { x: e.pos.x, y: e.pos.y, z: e.pos.z }, dist: frameDistance(e.radiusAU, 0.3) }, ang(v.yaw, v.pitch)); }
  if (name === 'orion') {
    const b = byId.betelgeuse, bp = b ? b.pos : { x: 0, y: 1e7, z: 0 }, bd = norm(bp);
    const side = norm(cross(bd, { x: 0, y: 0, z: 1 }));
    const v = norm({ x: side.x * 0.8 + 0.35 * 0 - bd.x * 0.45, y: side.y * 0.8 - bd.y * 0.45, z: 0.4 });
    const dist = Math.hypot(bp.x, bp.y, bp.z) * 1.15;
    return Object.assign({ pos: { x: bp.x * 0.45, y: bp.y * 0.45, z: bp.z * 0.45 }, dist: dist }, ang(Math.atan2(v.y, v.x) + 0.02 * t, Math.asin(clamp(v.z, -1, 1))));
  }
  if (name === 'galaxy') {
    const v = MW.obj.view, F = camFrame(v.yaw, v.pitch, 0);
    const th = TAU * Math.min(Y, 2.3e8) / 2.3e8 + 0.015 * t;  /* the camera circles the galaxy's axis, so the galaxy appears to turn */
    const d = rotateAbout(F.d, MW.ez, th), u = rotateAbout(F.up, MW.ez, th);
    const e = camAngles(d, u, v.yaw);
    return { pos: { x: MW.centre.x, y: MW.centre.y, z: MW.centre.z }, dist: objectViewDistance(MW.obj) * 1.1, yaw: e.yaw, pitch: e.pitch, roll: e.roll };
  }
  if (name === 'lg') {
    const c = MW.centre, m = byId.m31.pos, sep = Math.hypot(m.x - c.x, m.y - c.y, m.z - c.z);
    const v = byId.localgroup ? byId.localgroup.view : { yaw: 0, pitch: 0.5 };
    return Object.assign({ pos: { x: (c.x + m.x) / 2, y: (c.y + m.y) / 2, z: (c.z + m.z) / 2 }, dist: Math.max(sep * 1.25, 3.2e5 * LY_AU) }, ang(v.yaw + 0.02 * t, v.pitch));
  }
  if (name === 'cosmos') return Object.assign({ pos: { x: 0, y: 0, z: 0 }, dist: 6e8 * LY_AU }, ang(0.4 + 0.02 * t, 0.45));
  if (name === 'bh') { const b = byId.sgrA; return Object.assign({ pos: { x: b.pos.x, y: b.pos.y, z: b.pos.z }, dist: objectViewDistance(b) * 1.4 }, ang(0.3 + 0.05 * t, 0.25)); }
  return Object.assign({ pos: { x: 0, y: 0, z: 0 }, dist: 48 }, ang(-1.05, 0.42));
}
function deepShotName() {
  const L = deep.logY;
  if (L < 4.3) return 'solar';
  if (L < 6.6) return 'orion';
  if (L < 8.6) return 'galaxy';
  if (L < 9.15) return 'earth';
  if (L < 9.7) return 'lg';
  if (L < 10.3) return 'solarWide';
  if (L < 11.6) return 'cosmos';
  if (L < 30) return 'galaxy';
  return 'bh';
}
function deepCamera(dt) {
  const name = deepShotName();
  if (name !== deep.shot) {                                  /* start a smooth move from wherever the camera is */
    deep.shot = name; deep.blend = 0;
    deep.from = { pos: { x: cam.tx, y: cam.ty, z: cam.tz }, dist: cam.dist, yaw: cam.yaw, pitch: cam.pitch, roll: cam.roll };
  }
  const B = deepShot(name);
  let S = B;
  if (deep.blend < 1) {
    deep.blend = Math.min(1, deep.blend + dt / 2.2);
    const e = deep.blend < 0.5 ? 4 * deep.blend * deep.blend * deep.blend : 1 - Math.pow(-2 * deep.blend + 2, 3) / 2;
    const A = deep.from;
    const la = Math.log(A.dist), lb = Math.log(B.dist), d = Math.exp(lerp(la, lb, e));
    const sep = Math.hypot(B.pos.x - A.pos.x, B.pos.y - A.pos.y, B.pos.z - A.pos.z);
    let w = e;
    if (sep > 0) { const d0 = Math.max(A.dist, 0.3 * sep); if (B.dist > d0) w = clamp((Math.log(d) - Math.log(d0)) / (lb - Math.log(d0)), 0, 1); else if (A.dist > B.dist) w = clamp(1 - (Math.log(d) - Math.log(Math.max(B.dist, 0.3 * sep))) / (la - Math.log(Math.max(B.dist, 0.3 * sep)) || 1), 0, 1); }
    w = w * w * (3 - 2 * w);
    let dy = B.yaw - A.yaw; dy -= TAU * Math.round(dy / TAU);
    let dr = B.roll - A.roll; dr -= TAU * Math.round(dr / TAU);
    S = { pos: { x: lerp(A.pos.x, B.pos.x, w), y: lerp(A.pos.y, B.pos.y, w), z: lerp(A.pos.z, B.pos.z, w) }, dist: d,
          yaw: A.yaw + dy * e, pitch: lerp(A.pitch, B.pitch, e), roll: A.roll + dr * e };
  }
  cam.follow = null; cam.flight = null;
  cam.gx = cam.tx = S.pos.x; cam.gy = cam.ty = S.pos.y; cam.gz = cam.tz = S.pos.z;
  cam.dist = cam.distGoal = S.dist;
  cam.yaw = cam.yawGoal = S.yaw; cam.pitch = cam.pitchGoal = S.pitch; cam.roll = cam.rollGoal = S.roll;
  cam.pendYaw = cam.pendPitch = 0;
  cam.fov = cam.fovGoal = FOV_DEFAULT;
}
/* ---- what the renderer cannot show by itself ---- */
function deepOverlay() {
  const L = deep.logY;
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  /* Betelgeuse's supernova */
  const bet = byId.betelgeuse;
  if (bet && L >= 4.995 && L < 5.45) {
    const q = project(bet.pos);
    if (q) { const f = 1 - (L - 4.995) / 0.455; glow(q, 30 + 160 * f, 'rgb(220,235,255)', 0.9 * f); glow(q, 10 + 40 * f, 'rgb(255,255,255)', f); }
  }
  /* the Sun's planetary nebula */
  if (deep.Y > 7.72e9 && deep.Y < 8.8e9) {
    const q = project({ x: 0, y: 0, z: 0 });
    if (q) {
      const f = (deep.Y - 7.72e9) / 1.08e9, rAU = 1.5 + 9 * Math.sqrt(f), rp = (rAU / q.z) * focal, al = 0.55 * (1 - f);
      const g = ctx.createRadialGradient(q.x, q.y, rp * 0.6, q.x, q.y, rp);
      g.addColorStop(0, 'rgba(90,200,210,0)'); g.addColorStop(0.75, 'rgba(90,200,210,' + (al * 0.6).toFixed(3) + ')');
      g.addColorStop(0.93, 'rgba(255,120,150,' + al.toFixed(3) + ')'); g.addColorStop(1, 'rgba(255,120,150,0)');
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(q.x, q.y, rp, 0, TAU); ctx.fill();
    }
  }
  /* Milkomeda: the merged galaxy settles into a smooth elliptical glow */
  if (deep.merge > 0 && MW && deep.galaxyFade > 0.01) {
    const q = project(MW.centre);
    if (q) {
      const dC = Math.hypot(MW.centre.x - camPos.x, MW.centre.y - camPos.y, MW.centre.z - camPos.z);
      const rp = (MW.radius * 1.15 / dC) * focal, al = 0.42 * deep.merge * deep.galaxyFade;
      if (rp > 3) {
        const g = ctx.createRadialGradient(q.x, q.y, 0, q.x, q.y, rp);
        g.addColorStop(0, 'rgba(255,236,205,' + al.toFixed(3) + ')'); g.addColorStop(0.25, 'rgba(255,214,170,' + (al * 0.5).toFixed(3) + ')');
        g.addColorStop(1, 'rgba(255,200,160,0)');
        ctx.fillStyle = g; ctx.beginPath(); ctx.ellipse(q.x, q.y, rp, rp * 0.72, 0.4, 0, TAU); ctx.fill();
        if (deep.merge > 0.6 && deep.Y < 1e19) {
          ctx.globalCompositeOperation = 'source-over';
          ctx.font = '600 11px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'; ctx.textAlign = 'center';
          ctx.fillStyle = 'rgba(255,225,190,' + (0.8 * Math.min(1, deep.galaxyFade * 2)).toFixed(2) + ')';
          ctx.fillText('Milkomeda', q.x, q.y - rp * 0.78 - 8);
          ctx.globalCompositeOperation = 'lighter';
        }
      }
    }
  }
  /* Sagittarius A* evaporates in a last flash */
  const s = byId.sgrA;
  if (s && L >= 87.2 && L < 87.7) {
    const q = project(s.pos);
    if (q) { const f = 1 - (L - 87.2) / 0.5; glow(q, 40 + 260 * f, 'rgb(255,240,220)', f); }
  }
  ctx.restore();
}
/* ---- the readout ---- */
const DEEP_TAGS = { calc: 'Calculated', pred: 'Predicted', spec: 'Speculative' };
function deepHud() {
  const Y = deep.Y, t = fmtYears(Y);
  $('deepWhen').innerHTML = '+' + (t ? esc(t) : '10<sup>' + Math.floor(deep.logY) + '</sup> years');
  const yr = Math.round(deep.year0 + Y);
  $('deepSub').textContent = Y < 1e6 ? 'Year ' + (yr < 10000 ? String(yr) : yr.toLocaleString()) : 'from today';
  const r = deep.yrsPerSec;
  $('deepRate').textContent = !deep.playing ? 'Paused' : deep.hold > 0 ? 'Slowing for a moment' :
    'Time now passes at ' + (r < 1e15 && fmtYears(r) ? fmtYears(r) : '10^' + Math.floor(Math.log10(Math.max(r, 1))) + ' years') + ' per second';
  if (!deep.dragSlider) $('deepSlider').value = Math.round(deep.logY * 10);
  $('deepPlay').classList.toggle('on', deep.playing);
  $('deepPlayIcon').classList.toggle('hidden', deep.playing); $('deepPauseIcon').classList.toggle('hidden', !deep.playing);
  $('deepGuide').classList.toggle('hidden', !deep.takeover);
}
function deepShowChapter() {
  const c = DEEP_TIME[deep.chapter];
  const card = $('deepCard');
  if (!c) { card.classList.add('hidden'); return; }
  card.classList.remove('hidden');
  $('deepTag').textContent = DEEP_TAGS[c.tag]; $('deepTag').className = 'deep-tag ' + c.tag;
  $('deepTitle').textContent = c.title; $('deepText').textContent = c.text;
  card.classList.remove('pop'); void card.offsetWidth; card.classList.add('pop');
}
function deepChapterFor(L) { let k = -1; for (let i = 0; i < DEEP_TIME.length; i++) if (DEEP_TIME[i].log <= L + 1e-9) k = i; return k; }
function deepEnd(show) { $('deepEnd').classList.toggle('hidden', !show); }
/* ---- in and out ---- */
function deepStart() {
  if (deep.on) return;
  closeRoute(); closeCompare(); select(null);
  deepStore.panelOpen = !panel.classList.contains('hidden');
  deepStore.playing = playing; deepStore.simJd = simJd;
  deepStore.sun = { r: byId.sun.radiusAU, col: byId.sun.color };
  deepStore.earthColor = byId.earth.color; deepStore.earthColor2 = byId.earth.color2;
  deepStore.m31pos = { x: byId.m31.pos.x, y: byId.m31.pos.y, z: byId.m31.pos.z };
  panel.classList.add('hidden'); $('reopenBtn').classList.remove('show');
  $('timeBar').classList.add('hidden'); $('deepHud').classList.remove('hidden'); $('deepBar').classList.remove('hidden');
  deep.on = true; deep.jd0 = simJd; deep.year0 = dateFromJd(simJd).getUTCFullYear() + dateFromJd(simJd).getUTCMonth() / 12;
  deepRestart();
}
function deepRestart() {
  objects.forEach(function (o) { o.deepBlur = false; });
  deepSetLog(0); deep.chapter = -1; deep.hold = 0; deep.playing = true; deep.takeover = false; deep.shot = null; deep.clock = 0;
  deepEnd(false); $('deepCard').classList.add('hidden');
}
function deepStop() {
  if (!deep.on) return;
  deep.on = false;
  objects.forEach(function (o) { o.deepBlur = false; });
  byId.sun.radiusAU = deepStore.sun.r; byId.sun.color = deepStore.sun.col;
  byId.earth.dry = false; byId.earth.color = deepStore.earthColor; byId.earth.color2 = deepStore.earthColor2;
  byId.m31.pos = deepStore.m31pos;
  simJd = deepStore.simJd; playing = deepStore.playing; updatePlayIcon();
  $('timeBar').classList.remove('hidden'); $('deepHud').classList.add('hidden'); $('deepBar').classList.add('hidden'); deepEnd(false);
  cam.rollGoal = 0;
  if (deepStore.panelOpen) panel.classList.remove('hidden');
  goHome();
}
$('btnDeepTime').addEventListener('click', deepStart);
$('deepExit').addEventListener('click', deepStop);
$('deepEndExit').addEventListener('click', deepStop);
$('deepAgain').addEventListener('click', deepRestart);
$('deepPlay').addEventListener('click', function () {
  if (deep.logY >= 100.29) { deepRestart(); return; }
  deep.playing = !deep.playing; if (deep.playing) deepEnd(false);
});
$('deepGuide').addEventListener('click', function () { deep.takeover = false; deep.shot = null; });
const deepSlider = $('deepSlider');
deepSlider.addEventListener('input', function () {
  deep.dragSlider = true;
  const L = deepSlider.value / 10;
  objects.forEach(function (o) { o.deepBlur = false; });
  deepSetLog(L); deep.hold = 0; deep.chapter = deepChapterFor(L); deepShowChapter(); deepEnd(false);
});
deepSlider.addEventListener('change', function () { deep.dragSlider = false; });
/* chapter ticks along the slider */
(function () {
  const box = $('deepTicks');
  DEEP_TIME.forEach(function (c) {
    const t = document.createElement('span'); t.className = 'deep-tick'; t.style.left = (Math.min(c.log, 100) / 100.3 * 100) + '%'; t.title = c.title;
    box.appendChild(t);
  });
})();
/* dragging during deep time hands the camera to the viewer; time keeps running */
canvas.addEventListener('pointerdown', function () { if (deep.on) deep.takeover = true; });
canvas.addEventListener('wheel', function () { if (deep.on) deep.takeover = true; });

window.OU = { objects: objects, byId: byId, deep: deep, deepStart: deepStart, deepStop: deepStop, deepSetLog: deepSetLog, camFrame: camFrame, camAngles: camAngles, get view() { return { camPos: camPos, fwd: fwd, right: right, up: up }; }, openCompare: openCompare, bodyFrame: bodyFrame, surfVec: surfVec, discFrame: discFrame, spherePoly: spherePoly, deepList: deepList, MW: MW, STRUCT: STRUCT, cosmicAddress: cosmicAddress, galacticToEcl: galacticToEcl, cosmo: cosmo, COSMO: COSMO, lookAtDir: lookAtDir, openRoute: openRoute, closeRoute: closeRoute, route: route, fmtDuration: fmtDuration, cam: cam, opts: opts, flyTo: flyTo, select: select,
              /* render exactly one frame — useful when the tab is backgrounded
                 and requestAnimationFrame is throttled */
              step: function (dt, noDraw) { dt = dt || 0.016; if (deep.on) deepTick(dt); updatePositions(simJd); if (deep.on) deepWorld(); if (deep.on && !deep.takeover) deepCamera(dt); updateCamera(dt); if (!noDraw) { render(dt); if (deep.on) { deepOverlay(); deepHud(); } } },
              get jd() { return simJd; }, set jd(v) { simJd = v; } };

})();
