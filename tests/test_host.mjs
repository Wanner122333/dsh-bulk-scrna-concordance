import assert from 'node:assert/strict';
import { copyFile, mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { fileURLToPath, pathToFileURL } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

// Run a real child without a Python dependency. Node executes this .py entry
// file as JavaScript; the copied plugin still uses its normal script path.
const CHILD = String.raw`
const fs = require('node:fs');
const path = require('node:path');
const args = process.argv.slice(2);
const option = (name) => args[args.indexOf(name) + 1];
const outDir = option('--out-dir');
fs.writeFileSync(path.join(outDir, 'child-args.json'), JSON.stringify(args));
const mode = MODE;
if (mode === 'slow') Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, 200);
if (mode === 'fail') process.exit(1);
if (mode === 'no-write') process.exit(0);
const manifest = path.join(outDir, option('--out-prefix') + '_manifest.json');
fs.writeFileSync(manifest, '{}');
const result = {
  schema: 'dsh.bulk-scrna-concordance/1', status: 'ok',
  scrna: { n_cells: 2, n_genes_unique: 3 },
  pairs: { A: { n_cells: 2, n_genes_used: 3, pearson_r: 0.9,
    spearman_rho: 0.8, slope: 2, intercept: 1 } },
  warnings: [], outputs: { manifest_json: manifest }
};
fs.writeFileSync(option('--results-json'), JSON.stringify(result));
`;

async function fixture(t, mode = 'success') {
  const tempRoot = path.resolve(os.tmpdir());
  const dir = await mkdtemp(path.join(tempRoot, 'concordance-host-'));
  t.after(async () => {
    // Only remove the unique temporary directory created by this test.
    assert.equal(path.dirname(dir), tempRoot);
    assert.ok(path.basename(dir).startsWith('concordance-host-'));
    await rm(dir, { recursive: true, force: true });
  });
  await mkdir(path.join(dir, 'python'));
  await copyFile(path.join(ROOT, 'index.js'), path.join(dir, 'index.mjs'));
  await writeFile(path.join(dir, 'python', 'bulk_scrna_concordance.py'),
    CHILD.replace('MODE', JSON.stringify(mode)));
  const { apply } = await import(pathToFileURL(path.join(dir, 'index.mjs')).href);
  let tool;
  apply({
    effect(callback) { callback(); },
    tools: { register(value) { tool = value; } },
  }, { pythonPath: process.execPath, defaultTimeoutMs: 5000 });
  assert.ok(tool, 'the plugin registers its tool');
  const outDir = path.join(dir, 'out');
  await mkdir(outDir);
  const args = {
    sc_h5ad: path.join(dir, 'input.h5ad'),
    bulk_csv: path.join(dir, 'bulk.csv'),
    out_dir: outDir,
    out_prefix: 'test',
  };
  return {
    tool, args, outDir,
    childPath: path.join(dir, 'python', 'bulk_scrna_concordance.py'),
    execute: (overrides = {}) => tool.execute({ ...args, ...overrides }),
    childArgs: async () => JSON.parse(await readFile(path.join(outDir, 'child-args.json'), 'utf8')),
  };
}

async function seedOldResult(outDir) {
  await writeFile(path.join(outDir, 'test_results.json'), JSON.stringify({
    schema: 'dsh.bulk-scrna-concordance/1', status: 'ok',
    pairs: { OLD: { pearson_r: 1 } }, warnings: [], outputs: {},
  }));
}

test('a nonzero child cannot return a stale successful result', async (t) => {
  const f = await fixture(t, 'fail');
  await seedOldResult(f.outDir);
  const value = await f.execute();
  assert.equal(value.status, 'failed');
  assert.deepEqual(value.groups, []);
  assert.match(value.message, /code 1/);
});

test('a child that exits zero without writing cannot reuse old results', async (t) => {
  const f = await fixture(t, 'no-write');
  await seedOldResult(f.outDir);
  const value = await f.execute();
  assert.equal(value.status, 'failed');
  assert.deepEqual(value.groups, []);
});

test('fresh successful results expose their manifest path', async (t) => {
  const f = await fixture(t);
  const value = await f.execute();
  assert.equal(value.status, 'ok');
  assert.equal(value.groups[0].group, 'A');
  assert.equal(value.manifest_json, path.join(f.outDir, 'test_manifest.json'));
  assert.equal(await readFile(value.manifest_json, 'utf8'), '{}');
  const rendered = f.tool.output.render({}, value)[0].text;
  assert.ok(rendered.includes(value.manifest_json));
});

test('the selected counts layer is forwarded to Python', async (t) => {
  const f = await fixture(t);
  assert.equal((await f.execute({ sc_layer: 'counts' })).status, 'ok');
  const args = await f.childArgs();
  assert.equal(args[args.indexOf('--sc-layer') + 1], 'counts');
  assert.ok(!args.includes('--sc-use-raw'));
});

test('explicit raw selection is forwarded to Python', async (t) => {
  const f = await fixture(t);
  assert.equal((await f.execute({ sc_use_raw: true })).status, 'ok');
  const args = await f.childArgs();
  assert.ok(args.includes('--sc-use-raw'));
  assert.ok(!args.includes('--sc-layer'));
});

test('incompatible matrix selections fail before starting the child', async (t) => {
  const f = await fixture(t);
  await assert.rejects(f.execute({ sc_layer: 'counts', sc_use_raw: true }));
  await assert.rejects(readFile(path.join(f.outDir, 'child-args.json')), { code: 'ENOENT' });
});

for (const selection of [{ sc_layer: 'counts' }, { sc_use_raw: true }]) {
  const label = Object.keys(selection)[0];
  test(`CSV input rejects ${label} before starting the child`, async (t) => {
    const f = await fixture(t);
    await assert.rejects(f.execute({
      sc_h5ad: undefined, sc_matrix_csv: 'expression.csv', sc_meta_csv: 'cells.csv',
      ...selection,
    }));
    await assert.rejects(readFile(path.join(f.outDir, 'child-args.json')), { code: 'ENOENT' });
  });
}

test('one-group discriminability renders unavailable comparisons as NA', async (t) => {
  const f = await fixture(t);
  const value = await f.execute();
  value.discriminability = {
    matched_mean_r: 0.9, n_groups: 1, n_correct: 1, misassigned: [],
  };
  const rendered = f.tool.output.render({}, value)[0].text;
  assert.match(rendered, /mismatched=NA/);
  assert.match(rendered, /gap=NA/);
});

test('invalid output prefixes cannot start a child or remove an outside result', async (t) => {
  for (const prefix of ['../escape', 'sub/name', 'sub\\name', 'bad:name', 'bad\0name', '.', '..']) {
    await t.test(JSON.stringify(prefix), async (subtest) => {
      const f = await fixture(subtest);
      const outside = path.join(f.outDir, '..', 'escape_results.json');
      await writeFile(outside, 'keep this unrelated result');
      let rejection;
      try {
        await f.execute({ out_prefix: prefix });
      } catch (error) {
        rejection = error;
      }
      assert.equal(await readFile(outside, 'utf8'), 'keep this unrelated result');
      assert.ok(rejection instanceof Error, 'invalid prefix must be rejected');
      assert.match(rejection.message, /prefix/i);
      await assert.rejects(readFile(path.join(f.outDir, 'child-args.json')), { code: 'ENOENT' });
    });
  }
});

test('concurrent runs sharing an output target reject one caller without corrupting the result', async (t) => {
  const f = await fixture(t, 'slow');
  const values = await Promise.all([f.execute(), f.execute()]);
  const successes = values.filter((value) => value.status === 'ok');
  const failures = values.filter((value) => value.status === 'failed');
  assert.equal(successes.length, 1);
  assert.equal(failures.length, 1);
  assert.match(failures[0].message, /active|busy|already running/i);
  assert.equal(successes[0].groups[0].group, 'A');
  const result = JSON.parse(await readFile(path.join(f.outDir, 'test_results.json'), 'utf8'));
  assert.equal(result.pairs.A.pearson_r, 0.9);
});

test('a failed run releases its output target for a successful retry', async (t) => {
  const f = await fixture(t, 'fail');
  assert.equal((await f.execute()).status, 'failed');
  await writeFile(f.childPath, CHILD.replace('MODE', JSON.stringify('success')));
  const retry = await f.execute();
  assert.equal(retry.status, 'ok');
  assert.equal(retry.groups[0].group, 'A');
});

