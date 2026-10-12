// vswr.js
//
// Cascade of identical-topology stages, counted from the load towards the source:
//
//   Zload --[series line 1]--+--[series line 2]--+-- ... --[series line N]--+  entry point N
//                            |                   |                          |
//                         stub 1              stub 2                     stub N
//                     (shorted, shunt)    (shorted, shunt)           (shorted, shunt)
//
// Entry point k is where series line k meets stub k. The impedance seen at entry
// point k-1 becomes the load impedance of series line k (entry point 0 = Zload).
// The last entry point (N) is the entry point of the whole network; its reflection
// is evaluated against the system impedance Z0. All lines are lossless.
 
const C0 = 299792458; // m/s
 
// ---- small complex helpers: {re, im} ----
const cadd = (a, b) => ({ re: a.re + b.re, im: a.im + b.im });
const cabs = (a) => Math.hypot(a.re, a.im);
const cfinite = (a) => isFinite(a.re) && isFinite(a.im);
function cdiv(a, b) {
  const d = b.re * b.re + b.im * b.im;
  return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d };
}
const EPS = 1e-12;
const SHORT = { re: 0, im: 0 };
const OPEN = { re: Infinity, im: 0 };
 
/** Electrical length (rad) of a line of physical length `len` (m) at `freqMHz`. */
function theta(freqMHz, len, vf) {
  return 2 * Math.PI * freqMHz * 1e6 * len / (C0 * vf);
}
 
/**
 * Series line (Z0b, theta) terminated in ZL, seen from the entry point:
 *   Zin = Z0b (ZL cos + j Z0b sin) / (Z0b cos + j ZL sin)
 * Returns { Y, Zin } where Y = 1/Zin (computed so that short/open cases stay safe).
 */
function seriesLine(ZL, Z0b, th) {
  const c = Math.cos(th), s = Math.sin(th);
 
  if (!cfinite(ZL)) {                                   // open-circuited load: Zin = -j Z0b cot(theta)
    if (Math.abs(s) < EPS) return { Y: { re: 0, im: 0 }, Zin: OPEN };
    if (Math.abs(c) < EPS) return { Y: { re: Infinity, im: 0 }, Zin: SHORT };
    return { Y: { re: 0, im: s / (Z0b * c) }, Zin: { re: 0, im: -Z0b * c / s } };
  }
 
  const num = { re: Z0b * ZL.re * c, im: Z0b * (ZL.im * c + Z0b * s) };   // Z0b (ZL c + j Z0b s)
  const den = { re: Z0b * c - ZL.im * s, im: ZL.re * s };                  // Z0b c + j ZL s
  if (cabs(num) < EPS) return { Y: { re: Infinity, im: 0 }, Zin: SHORT };
  if (cabs(den) < EPS) return { Y: { re: 0, im: 0 }, Zin: OPEN };
  return { Y: cdiv(den, num), Zin: cdiv(num, den) };
}
 
/** Admittance of a shorted stub: Y = 1/(j Z0a tan th) = -j cos / (Z0a sin). */
function shortedStubY(Z0a, th) {
  const c = Math.cos(th), s = Math.sin(th);
  if (Math.abs(s) < EPS) return { re: Infinity, im: 0 };   // theta = n*pi: stub looks like a short
  return { re: 0, im: -c / (Z0a * s) };
}
 
/** Reflection coefficient of impedance Z (complex) against real Z0. */
function reflectionCoefficient(Z, Z0) {
  if (!cfinite(Z)) return { re: 1, im: 0, mag: 1 };
  const g = cdiv({ re: Z.re - Z0, im: Z.im }, { re: Z.re + Z0, im: Z.im });
  return { re: g.re, im: g.im, mag: cabs(g) };
}
 
/** VSWR = (1 + |G|) / (1 - |G|), Infinity when |G| >= 1. */
function vswrFromGamma(mag) {
  if (!isFinite(mag) || mag >= 1 - 1e-12) return Infinity;
  return (1 + mag) / (1 - mag);
}
 
/**
 * One stage at one frequency: series line terminated in `Zload`, then the shunt
 * stub at the stage's entry point.
 * stage: { line: {Z0, len, vf}, stub: {enabled, Z0, len, vf} }
 * Returns { Zline, Zentry } (Zline = impedance after the series line, before the stub).
 */
function stageImpedance(freq, Zload, stage) {
  const ln = seriesLine(Zload, stage.line.Z0, theta(freq, stage.line.len, stage.line.vf));
  let Y = ln.Y;
 
  if (stage.stub.enabled) {
    const ys = shortedStubY(stage.stub.Z0, theta(freq, stage.stub.len, stage.stub.vf));
    if (!cfinite(ys)) return { Zline: ln.Zin, Zentry: SHORT };   // stub shorts the entry point
    Y = cadd(Y, ys);
  }
 
  if (!cfinite(Y)) return { Zline: ln.Zin, Zentry: SHORT };
  if (cabs(Y) < EPS) return { Zline: ln.Zin, Zentry: OPEN };
  return { Zline: ln.Zin, Zentry: cdiv({ re: 1, im: 0 }, Y) };
}
 
/**
 * Impedance at the entry point of the whole cascade for one frequency.
 * cfg: { Z0, stages: [stage, stage, ...] }  (stages[0] is next to the load)
 * Returns { entries: [Zentry of stage 1, ..., Zentry of stage N], Zentry }.
 */
function entryImpedance(freq, Zload, cfg) {
  let Z = Zload;
  const entries = [];
  for (const stage of cfg.stages) {
    Z = stageImpedance(freq, Z, stage).Zentry;   // previous entry point becomes the next Zload
    entries.push(Z);
  }
  return { entries, Zentry: Z };
}
 
/**
 * points: [{freq (MHz), R, X}]  -> Zload = R + jX
 * returns { results, max }
 */
function analyze(points, cfg) {
  const results = points.map((p, i) => {
    const Zload = { re: p.R, im: p.X };
    const e = entryImpedance(p.freq, Zload, cfg);
    const g = reflectionCoefficient(e.Zentry, cfg.Z0);
    return {
      index: i,
      freq: p.freq,
      Zload: Zload,
      entries: e.entries,      // entry-point impedance after each stage
      Zentry: e.Zentry,        // final entry point of the whole network
      gammaMag: g.mag,
      gammaAngleDeg: Math.atan2(g.im, g.re) * 180 / Math.PI,
      returnLossDb: g.mag === 0 ? Infinity : -20 * Math.log10(g.mag),
      vswr: vswrFromGamma(g.mag)
    };
  });
  let max = null;
  for (const r of results) if (max === null || r.vswr > max.vswr) max = r;
  return { results, max };
}
 
if (typeof module !== "undefined") {
  module.exports = { analyze, entryImpedance, stageImpedance, reflectionCoefficient, vswrFromGamma };
}
 
// ================= UI (browser only) =================
if (typeof document !== "undefined") {
  const $ = (id) => document.getElementById(id);
  const tbody = $("rows"), resultsBody = $("results"), resultsHead = $("resHead");
  const summary = $("summary"), errorBox = $("error"), stagesBox = $("stages");
 
  // ---------- Z_load table: each click on "Add frequency" appends a blank row ----------
  function addRow(p) {
    p = p || { freq: "", R: "", X: "" };
    const tr = document.createElement("tr");
    tr.innerHTML =
      '<td><input type="number" step="any" class="f" value="' + p.freq + '"></td>' +
      '<td><input type="number" step="any" class="r" value="' + p.R + '"></td>' +
      '<td><input type="number" step="any" class="x" value="' + p.X + '"></td>' +
      '<td><button type="button" class="del" title="Remove row">&times;</button></td>';
    tr.querySelector(".del").addEventListener("click", () => { tr.remove(); calculate(); });
    tr.querySelectorAll("input").forEach((i) => i.addEventListener("input", calculate));
    tbody.appendChild(tr);
    tr.querySelector(".f").focus();
  }
 
  // ---------- stages: each "Add stage" appends a series line + a shunt stub ----------
  function stageNumberLabels() {
    const stages = stagesBox.querySelectorAll(".stage");
    stages.forEach((st, i) => {
      const n = i + 1;
      const prev = n === 1 ? "Z<sub>load</sub>" : "entry point " + (n - 1);
      st.querySelector(".stageTitle").innerHTML =
        "Stage " + n + ": series line " + n + " (between " + prev + " and entry point " + n +
        ") + shorted stub " + n + " at entry point " + n;
      st.querySelector(".delStage").style.display = stages.length > 1 ? "" : "none";
    });
  }
 
  function addStage() {
    const st = document.createElement("div");
    st.className = "stage";
    st.innerHTML =
      '<div class="stageHead"><strong class="stageTitle"></strong>' +
      '<button type="button" class="delStage" title="Remove this stage">&times; Remove stage</button></div>' +
      '<div class="cfg">' +
      '<fieldset><legend>Series line</legend>' +
      '<label>Characteristic impedance (&Omega;) <input type="number" class="lineZ0" value="50" step="any" min="0"></label>' +
      '<label>Physical length (m) <input type="number" class="lineLen" value="0.15" step="any" min="0"></label>' +
      '<label>Velocity factor <input type="number" class="lineVf" value="0.66" step="any" min="0" max="1"></label>' +
      '</fieldset>' +
      '<fieldset><legend>Shunt stub, short-circuited at far end</legend>' +
      '<label>Include stub <input type="checkbox" class="stubOn" checked></label>' +
      '<label>Characteristic impedance (&Omega;) <input type="number" class="stubZ0" value="50" step="any" min="0"></label>' +
      '<label>Physical length (m) <input type="number" class="stubLen" value="0.10" step="any" min="0"></label>' +
      '<label>Velocity factor <input type="number" class="stubVf" value="0.66" step="any" min="0" max="1"></label>' +
      '</fieldset></div>';
    st.querySelector(".delStage").addEventListener("click", () => {
      st.remove(); stageNumberLabels(); calculate();
    });
    st.querySelectorAll("input").forEach((i) => {
      i.addEventListener("input", calculate);
      i.addEventListener("change", calculate);
    });
    stagesBox.appendChild(st);
    stageNumberLabels();
    calculate();
  }
 
  const val = (el, sel) => parseFloat(el.querySelector(sel).value);
 
  function readStages() {
    return Array.prototype.map.call(stagesBox.querySelectorAll(".stage"), (st) => ({
      line: { Z0: val(st, ".lineZ0"), len: val(st, ".lineLen"), vf: val(st, ".lineVf") },
      stub: {
        enabled: st.querySelector(".stubOn").checked,
        Z0: val(st, ".stubZ0"), len: val(st, ".stubLen"), vf: val(st, ".stubVf")
      }
    }));
  }
 
  function fmt(v, d) {
    if (v === Infinity) return "∞";
    return Number(v).toFixed(d === undefined ? 3 : d);
  }
  function fmtZ(z, d) {
    d = d === undefined ? 2 : d;
    if (!isFinite(z.re) || !isFinite(z.im)) return "∞";
    return fmt(z.re, d) + (z.im < 0 ? " - j" : " + j") + fmt(Math.abs(z.im), d);
  }
 
  function calculate() {
    errorBox.textContent = "";
    const cfg = { Z0: parseFloat($("z0").value), stages: readStages() };
 
    const bad = [];
    if (!(cfg.Z0 > 0)) bad.push("Z0");
    const lineOk = (l) => l.Z0 > 0 && l.len >= 0 && l.vf > 0 && l.vf <= 1;
    cfg.stages.forEach((s, i) => {
      if (!lineOk(s.line)) bad.push("stage " + (i + 1) + " series line (Z0 > 0, length >= 0, 0 < VF <= 1)");
      if (s.stub.enabled && !lineOk(s.stub)) bad.push("stage " + (i + 1) + " stub (Z0 > 0, length >= 0, 0 < VF <= 1)");
    });
    if (bad.length) { errorBox.textContent = "Check inputs: " + bad.join("; "); return; }
 
    const points = [];
    tbody.querySelectorAll("tr").forEach((tr) => {
      const f = parseFloat(tr.querySelector(".f").value);
      const R = parseFloat(tr.querySelector(".r").value);
      const X = parseFloat(tr.querySelector(".x").value);
      if (!isNaN(f) && f > 0 && !isNaN(R) && !isNaN(X)) points.push({ freq: f, R, X });
    });
 
    // header: one Z_entry column per stage
    let head = "<tr><th>f (MHz)</th><th>Z<sub>load</sub> (&Omega;)</th>";
    cfg.stages.forEach((s, i) => { head += "<th>Z<sub>entry " + (i + 1) + "</sub> (&Omega;)</th>"; });
    head += "<th>|&Gamma;|</th><th>Return loss (dB)</th><th>VSWR</th></tr>";
    resultsHead.innerHTML = head;
    resultsBody.innerHTML = "";
 
    if (!points.length) {
      summary.textContent = "Click \"Add frequency\", then enter the frequency and Zload (Re, Im) in the new row.";
      return;
    }
 
    const out = analyze(points, cfg);
    out.results.forEach((r) => {
      const tr = document.createElement("tr");
      if (r.index === out.max.index) tr.className = "max";
      let html = "<td>" + r.freq + "</td><td>" + fmtZ(r.Zload) + "</td>";
      r.entries.forEach((z) => { html += "<td>" + fmtZ(z) + "</td>"; });
      html += "<td>" + fmt(r.gammaMag, 4) + "</td><td>" + fmt(r.returnLossDb, 2) + "</td><td>" + fmt(r.vswr, 3) + "</td>";
      tr.innerHTML = html;
      resultsBody.appendChild(tr);
    });
    summary.innerHTML = "Maximum VSWR = <strong>" + fmt(out.max.vswr, 3) +
      "</strong> at <strong>" + out.max.freq + " MHz</strong> (Z<sub>entry " + cfg.stages.length + "</sub> = " +
      fmtZ(out.max.Zentry) + " &Omega;)";
  }
 
  $("add").addEventListener("click", () => addRow());
  $("addStage").addEventListener("click", () => addStage());
  $("z0").addEventListener("input", calculate);
 
  addStage();   // stage 1 (series line + stub) is always present
}
 