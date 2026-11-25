import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { Workbook } from "@oai/artifact-tool";

const directory = process.argv[2];
const sourceKeys = ["case", "U_m_s", "alpha_deg", "sweep_deg", "eta"];
const key = (row) => JSON.stringify(sourceKeys.map((name) => Number(row[name])));

async function readCsv(name) {
  const workbook = await Workbook.fromCSV(await fs.readFile(path.join(directory, name), "utf8"), {
    sheetName: "Source",
  });
  const [header, ...rows] = workbook.worksheets.getItem("Source").getUsedRange(true).values;
  return rows.map((row) => Object.fromEntries(header.map((name, column) => [name, row[column]])));
}

function indexRows(rows) {
  const indexed = new Map(rows.map((row) => [key(row), row]));
  assert.equal(indexed.size, rows.length);
  return indexed;
}

const [forces, perturbations, moments] = await Promise.all([
  readCsv("spanwise_cn_alpha.csv"),
  readCsv("spanwise_cn_alpha_perturbations.csv"),
  readCsv("spanwise_cm.csv"),
]);
const fixedShapeByStation = indexRows(perturbations.filter((row) => row.mode === "fixed_shape"));
const momentsByStation = indexRows(moments);
assert.equal(indexRows(forces).size, 4);
assert.equal(fixedShapeByStation.size, forces.length);
assert.equal(momentsByStation.size, forces.length);

const header = [
  "U_m_s", "alpha_deg", "sweep_deg", "eta", "alpha_step_deg", "alpha_local_deg", "cn",
  "cn_alpha_secant_per_rad", "cn_alpha_fixed_shape_local_per_rad", "fixed_shape_delta_alpha_local_rad",
  "cm_quarter_chord", "cm_spar",
];
const cleanRows = forces.map((force) => {
  const fixedShape = fixedShapeByStation.get(key(force));
  const moment = momentsByStation.get(key(force));
  assert(fixedShape && moment, "Each span station needs matching fixed-shape and moment results");
  assert.equal(String(force.central_converged), "true");
  assert.equal(String(force.cn_alpha_local_valid), "true");
  assert.equal(String(fixedShape.valid), "true");
  assert.equal(String(moment.central_converged), "true");
  assert.equal(Number(fixedShape.alpha_step_deg), Number(force.alpha_step_deg));
  assert(Math.abs(Number(fixedShape.delta_cn) / Number(fixedShape.delta_alpha_local_rad) -
    Number(fixedShape.cn_alpha_local_perturbation_per_rad)) < 1e-12);
  const row = [
    force.U_m_s, force.alpha_deg, force.sweep_deg, force.eta, force.alpha_step_deg, force.alpha_local_deg,
    force.cn, force.cn_alpha_local_per_rad, fixedShape.cn_alpha_local_perturbation_per_rad,
    fixedShape.delta_alpha_local_rad, moment.cm_quarter_chord, moment.cm_spar,
  ].map(String);
  assert(row.every((value) => value !== "" && Number.isFinite(Number(value))));
  return row;
});

const output = path.join(directory, "spanwise_cn_alpha_clean.csv");
const content = [header, ...cleanRows].map((row) => row.join(",")).join("\n") + "\n";
const check = await Workbook.fromCSV(content, { sheetName: "Spanwise coefficients" });
assert.deepEqual(check.worksheets.getItemAt(0).getUsedRange(true).values, [header, ...cleanRows]);
await fs.writeFile(output, content, { encoding: "utf8", flag: "wx" });
assert.equal(await fs.readFile(output, "utf8"), content);
console.log(`Saved and verified ${cleanRows.length} stations and ${header.length} fields in ${output}`);
