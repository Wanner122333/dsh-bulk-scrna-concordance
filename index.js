/**
 * Bulk–scRNA concordance analysis — DSH Host plugin.
 *
 * Registers one agent-callable tool, `bulk_scrna_concordance`, that runs a
 * local Python analysis script and returns a structured, size-bounded result.
 *
 * Two decisions worth knowing about:
 *
 * 1. The computation stays in Python. The original tool is Python (anndata /
 *    h5py / scipy); reading `.h5ad` in Node would mean a second, differently
 *    behaved implementation of the same statistics. The plugin is the thin,
 *    contract-owning layer; the script owns the numerics.
 *
 * 2. Child stdout goes to a log FILE, not a pipe. The Harness file sandbox
 *    denies piped stdio to child processes (`spawnSync ... EPERM`), while
 *    `stdio: 'inherit'` and fd redirection both work. Redirecting to a file
 *    keeps the output readable for diagnostics without any pipe.
 *
 * 3. This plugin imports NOTHING outside `node:`. A Host-only bundle needs no
 *    dependencies, and declaring one would make the install require the npm
 *    registry. The tool definition is therefore a plain registry-ready object
 *    with JSON Schema parameters, plus the local argument check that
 *    `defineTool` would otherwise supply.
 *
 * @module @local/dsh-bulk-scrna-concordance
 */

import { spawn } from 'node:child_process';
import { closeSync, openSync } from 'node:fs';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SCRIPT = path.join(HERE, 'python', 'bulk_scrna_concordance.py');
const RESULTS_SCHEMA = 'dsh.bulk-scrna-concordance/1';

/**
 * Runtime configuration keys this plugin understands.
 *
 * Deliberately NOT exported as `Config`. Cordis resolves an exported `Config`
 * through its schema library (`resolveConfig` calls `Config.validate`), and that
 * library cannot be resolved from this plugin's directory when the bundle is
 * installed by path -- importing it would make the plugin fail to load with
 * "Cannot read properties of undefined (reading 'validate')". `Config` is an
 * optional export, so it is omitted and the row's `config` values are read
 * defensively in `apply()` instead.
 */
const CONFIG_KEYS = ['pythonPath', 'defaultTimeoutMs', 'maxLogLines'];

/**
 * Read the plugin row's config with defaults.
 *
 * @param config - the `config` object from the bundle's patch row.
 * @returns resolved settings, replacing malformed values with their defaults.
 */
function resolveConfig(config) {
	const raw = config !== null && typeof config === 'object' ? config : {};
	const unknown = Object.keys(raw).filter((key) => !CONFIG_KEYS.includes(key));
	if (unknown.length > 0) {
		// A warning, not a failure: an unknown key must never stop activation.
		console.warn(`[bulk-scrna-concordance] ignoring unknown config key(s): ${unknown.join(', ')}`);
	}
	const pythonPath = typeof raw.pythonPath === 'string' && raw.pythonPath.trim() !== ''
		? raw.pythonPath.trim()
		: 'python';
	const defaultTimeoutMs = Number.isFinite(raw.defaultTimeoutMs) && raw.defaultTimeoutMs > 0
		? raw.defaultTimeoutMs
		: 1_800_000;
	const maxLogLines = Number.isInteger(raw.maxLogLines) && raw.maxLogLines > 0
		? raw.maxLogLines
		: 40;
	return { pythonPath, defaultTimeoutMs, maxLogLines };
}

export const name = 'bulk-scrna-concordance';
export const inject = ['tools'];

/**
 * Validate arguments against the declared parameter schema.
 *
 * This is the check `defineTool` would perform. Only the constructs this
 * plugin's schema actually uses are implemented: type, required, enum, and
 * per-item object properties.
 *
 * @returns An array of human-readable violations; empty means valid.
 */
function validateArgs(schema, args) {
	const violations = [];
	if (args === null || typeof args !== 'object' || Array.isArray(args)) {
		return ['arguments must be an object'];
	}
	for (const key of schema.required ?? []) {
		if (args[key] === undefined || args[key] === null) {
			violations.push(`missing required argument: ${key}`);
		}
	}
	for (const [key, spec] of Object.entries(schema.properties ?? {})) {
		const value = args[key];
		if (value === undefined || value === null) continue;
		const expected = spec.type;
		const actual = Array.isArray(value) ? 'array' : typeof value;
		const ok = expected === 'integer'
			? Number.isInteger(value)
			: expected === 'number'
				? typeof value === 'number' && Number.isFinite(value)
				: actual === expected;
		if (!ok) {
			violations.push(`${key} must be ${expected}, got ${actual}`);
			continue;
		}
		if (Array.isArray(spec.enum) && !spec.enum.includes(value)) {
			violations.push(`${key} must be one of ${spec.enum.join(', ')}`);
		}
	}
	return violations;
}

function requireNonEmpty(value, label) {
	if (typeof value !== 'string' || value.trim() === '') {
		throw new Error(`${label} must be a non-empty string`);
	}
	return value;
}

/**
 * Resolve a caller-supplied path against the session working directory.
 * Absolute paths pass through unchanged.
 */
function resolvePath(value, cwd) {
	const expanded = value.startsWith('~') ? path.join(os.homedir(), value.slice(1)) : value;
	return path.resolve(cwd ?? process.cwd(), expanded);
}

/**
 * Run the analysis script and collect its outcome.
 *
 * Returns `{ code, signal, logPath }`. Never rejects for a non-zero exit: the
 * caller distinguishes "ran and failed" from "could not start".
 */
function runScript(argv, { pythonPath, timeoutMs, logPath, cwd }) {
	return new Promise((resolve, reject) => {
		let settled = false;
		let timer;

		// stdout/stderr go to a raw file descriptor. `spawn` rejects a
		// WriteStream whose fd is not yet assigned, and a pipe is unavailable
		// under the file sandbox. The child inherits a duplicate of this fd, so
		// the parent can close its own copy immediately after spawning.
		let logFd;
		try {
			logFd = openSync(logPath, 'w');
		} catch (error) {
			reject(error);
			return;
		}

		let child;
		try {
			child = spawn(pythonPath, argv, {
				cwd,
				stdio: ['ignore', logFd, logFd],
				windowsHide: true,
			});
		} catch (error) {
			closeSync(logFd);
			reject(error);
			return;
		}
		closeSync(logFd);

		const finish = (payload) => {
			if (settled) return;
			settled = true;
			if (timer !== undefined) clearTimeout(timer);
			resolve({ ...payload, logPath });
		};

		child.on('error', (error) => {
			if (settled) return;
			settled = true;
			if (timer !== undefined) clearTimeout(timer);
			reject(error);
		});

		child.on('close', (code, signal) => finish({ code, signal: signal ?? null }));

		if (Number.isFinite(timeoutMs) && timeoutMs > 0) {
			timer = setTimeout(() => {
				try {
					child.kill();
				} catch {
					// already gone
				}
			}, timeoutMs);
		}
	});
}

/** Read the tail of a text file without loading all of it when it is large. */
async function readTail(file, maxLines) {
	try {
		const raw = await fs.readFile(file, 'utf8');
		const lines = raw.split(/\r?\n/).filter((line) => line.length > 0);
		return { lines: lines.slice(-maxLines), total: lines.length };
	} catch {
		return { lines: [], total: 0 };
	}
}

export function apply(ctx, config = {}) {
	const { pythonPath, defaultTimeoutMs, maxLogLines } = resolveConfig(config);

	const parameters = {
		type: 'object',
		additionalProperties: false,
		required: ['bulk_csv', 'out_dir'],
		properties: {
			sc_h5ad: {
				type: 'string',
				description: 'Path to the single-cell AnnData .h5ad file. Use this or sc_matrix_csv.',
			},
			sc_matrix_csv: {
				type: 'string',
				description: 'Path to a single-cell expression matrix (genes as rows, cells as columns) as an alternative to sc_h5ad.',
			},
			sc_meta_csv: {
				type: 'string',
				description: 'Cell metadata for sc_matrix_csv; first column is the cell ID and must include the group column.',
			},
			bulk_csv: {
				type: 'string',
				description: 'Bulk expression matrix (CSV/TSV/TXT); first column is the gene ID, remaining columns are samples.',
			},
			out_dir: {
				type: 'string',
				description: 'Directory for the result files. Created if it does not exist.',
			},
			out_prefix: {
				type: 'string',
				description: 'Filename prefix for the outputs. Defaults to bulk_scrna_concordance.',
			},
			group_col: {
				type: 'string',
				description: 'Metadata column holding the group label (L1/L2/L3 style). Defaults to Group.',
			},
			gene_id_key: {
				type: 'string',
				description: "Which var column holds gene IDs in the .h5ad, or 'auto' or 'var_names'. Defaults to auto.",
			},
			strip_version: {
				type: 'boolean',
				description: 'Strip the suffix after the first "." from gene IDs on both sides. Defaults to false. Turn on only when IDs genuinely carry a version suffix.',
			},
			gene_id_case: {
				type: 'string',
				enum: ['keep', 'upper', 'lower'],
				description: 'Case-normalize gene IDs on both sides before matching. Defaults to keep.',
			},
			discriminability: {
				type: 'boolean',
				description: 'Also compute the full pseudo-bulk x bulk correlation matrix with matched-vs-mismatched summary. Recommended.',
			},
			normalize: {
				type: 'string',
				enum: ['cpm', 'tpm', 'none'],
				description: "Library-size scaling before log2(x+1). 'cpm' is the default and the historical behaviour; 'tpm' also divides by gene length and requires sc_lengths_col and bulk_lengths_col; 'none' assumes the inputs are already on a common scale.",
			},
			sc_lengths_col: {
				type: 'string',
				description: 'AnnData var column holding gene length, for normalize=tpm.',
			},
			bulk_lengths_col: {
				type: 'string',
				description: 'Bulk column holding gene length, for normalize=tpm. The column is excluded from the expression matrix.',
			},
			no_manifest: {
				type: 'boolean',
				description: 'Skip manifest.json. By default every run writes one: input and output SHA-256 hashes, parameters, seed, and a pinned environment lock, so a published figure can be traced back to exact inputs.',
			},
			no_figure: {
				type: 'boolean',
				description: 'Skip the PNG/PDF/SVG figure and produce tables only.',
			},
			timeout_ms: {
				type: 'integer',
				description: 'Kill the analysis after this many milliseconds. Defaults to the plugin setting.',
			},
		},
	};

	// Output schema: `required` is a plain JSON Schema list here, because this
	// object is registry-ready rather than DSL input for defineTool().
	const outputSchema = {
		type: 'object',
		additionalProperties: false,
		required: ['status', 'groups'],
		properties: {
			status: { type: 'string', description: 'ok, failed, or killed.' },
			message: { type: 'string', description: 'Failure detail when status is not ok.' },
			n_cells: { type: 'integer', description: 'Cells in the single-cell input.' },
			n_genes_unique: { type: 'integer', description: 'Genes after collapsing duplicate IDs.' },
			groups: {
				type: 'array',
				description: 'Per-group concordance statistics.',
				items: {
					type: 'object',
					additionalProperties: false,
					required: ['group', 'pearson_r'],
					properties: {
						group: { type: 'string' },
						n_cells: { type: 'integer' },
						n_genes_used: { type: 'integer' },
						pearson_r: { type: 'number' },
						spearman_rho: { type: 'number' },
						slope: { type: 'number' },
						intercept: { type: 'number' },
					},
				},
			},
			discriminability: {
				type: 'object',
				additionalProperties: false,
				description: 'Present only when discriminability was requested.',
				properties: {
					matched_mean_r: { type: 'number' },
					mismatched_mean_r: { type: 'number' },
					gap: { type: 'number' },
					n_groups: { type: 'integer' },
					n_correct: { type: 'integer' },
					misassigned: {
						type: 'array',
						items: { type: 'string' },
						description: 'Pseudo-bulk groups whose best-matching bulk group is not their own.',
					},
				},
			},
					warnings: { type: 'array', items: { type: 'string' } },
					results_json: { type: 'string' },
					manifest_json: { type: 'string' },
					stats_csv: { type: 'string' },
					figure_png: { type: 'string' },
					figure_pdf: { type: 'string' },
					log_tail: { type: 'array', items: { type: 'string' } },
				},
	};

	const tool = {
		name: 'bulk_scrna_concordance',
		description: [
			'Compare bulk RNA-seq expression against scRNA-seq pseudo-bulk expression and report concordance.',
			'',
			'For each group it sums single-cell counts into a pseudo-bulk profile, CPM-normalizes both sides,',
			'applies log2(x+1), and computes Pearson/Spearman correlation and a linear fit over the shared genes.',
			'',
			'Set `discriminability: true` to also correlate every pseudo-bulk group against every bulk group.',
			'Do that whenever the correlation is meant as evidence: a single high r cannot distinguish a',
			'matching tissue from a mismatched one, because any two transcriptomes from the same experiment',
			'are strongly correlated. The discriminability block reports the matched-vs-mismatched gap and',
			'whether each pseudo-bulk group actually best-matches its own bulk group.',
			'',
			'Relative paths resolve against the session working directory. Returns a compact summary plus the',
			'paths of the full results JSON, statistics CSV, and figure.',
		].join('\n'),
		parameters,
		output: {
			schema: outputSchema,
			render: (_args, value) => {
				const lines = [];
				if (value.status !== 'ok') {
					lines.push(`bulk–scRNA concordance: ${value.status}`);
					if (value.message !== undefined) lines.push(value.message);
				} else {
					lines.push('bulk–scRNA concordance');
					for (const g of value.groups ?? []) {
						lines.push(
							`  ${g.group}: Pearson r=${g.pearson_r.toFixed(3)}`
							+ ` Spearman rho=${(g.spearman_rho ?? 0).toFixed(3)}`
							+ ` slope=${(g.slope ?? 0).toFixed(3)}`
							+ ` n_genes=${g.n_genes_used ?? 0}`
							+ (g.n_cells === undefined ? '' : ` n_cells=${g.n_cells}`),
						);
					}
					const d = value.discriminability;
					if (d !== undefined) {
						lines.push(
							`  discriminability: matched=${(d.matched_mean_r ?? 0).toFixed(4)}`
							+ ` mismatched=${(d.mismatched_mean_r ?? 0).toFixed(4)}`
							+ ` gap=${(d.gap ?? 0).toFixed(4)}`
							+ ` correct=${d.n_correct ?? 0}/${d.n_groups ?? 0}`,
						);
						if ((d.misassigned ?? []).length > 0) {
							lines.push(`  misassigned: ${d.misassigned.join(', ')}`);
						}
					}
				}
				for (const w of value.warnings ?? []) lines.push(`  warning: ${w}`);
				if (value.manifest_json !== undefined) lines.push(`  manifest: ${value.manifest_json}`);
				if (value.results_json !== undefined) lines.push(`  results: ${value.results_json}`);
				if (value.stats_csv !== undefined) lines.push(`  stats:   ${value.stats_csv}`);
				// The figure cannot be carried inline: Harness production adapters
				// declare text-only assistant output, so only user content may hold
				// an image block. Emit a Markdown image reference so the path is
				// clickable in the transcript, and tell the caller that surfacing
				// the figure needs the separate deliverable-card mechanism.
				if (value.figure_png !== undefined) {
					lines.push(`  figure:  ${value.figure_png}`);
					lines.push(`  figure (markdown): ![bulk–scRNA concordance](${value.figure_png})`);
					if (value.figure_pdf !== undefined) lines.push(`  figure (pdf): ${value.figure_pdf}`);
					lines.push('  note: the figure is a file on disk; to show it to the user as a card, call present with its path.');
				}
				if (value.status !== 'ok' && (value.log_tail ?? []).length > 0) {
					lines.push('  log tail:');
					for (const line of value.log_tail) lines.push(`    ${line}`);
				}
				return [{ type: 'text', text: lines.join('\n') }];
			},
		},
		async execute(args, exec) {
			const violations = validateArgs(parameters, args);
			if (violations.length > 0) throw new Error(violations.join('; '));

			const cwd = exec?.agent?.session?.cwd ?? process.cwd();
			const outDir = resolvePath(requireNonEmpty(args.out_dir, 'out_dir'), cwd);
			const bulkCsv = resolvePath(requireNonEmpty(args.bulk_csv, 'bulk_csv'), cwd);

			const hasH5ad = typeof args.sc_h5ad === 'string' && args.sc_h5ad.trim() !== '';
			const hasMatrix = typeof args.sc_matrix_csv === 'string' && args.sc_matrix_csv.trim() !== '';
			if (hasH5ad === hasMatrix) {
				throw new Error('provide exactly one of sc_h5ad or sc_matrix_csv');
			}
			if (hasMatrix && (typeof args.sc_meta_csv !== 'string' || args.sc_meta_csv.trim() === '')) {
				throw new Error('sc_matrix_csv requires sc_meta_csv');
			}

			const prefix = typeof args.out_prefix === 'string' && args.out_prefix.trim() !== ''
				? args.out_prefix.trim()
				: 'bulk_scrna_concordance';
			const resultsJson = path.join(outDir, `${prefix}_results.json`);
			const logPath = path.join(outDir, `${prefix}_run.log`);

			await fs.mkdir(outDir, { recursive: true });

			try {
				await fs.access(SCRIPT);
			} catch {
				throw new Error(`analysis script missing from the plugin package: ${SCRIPT}`);
			}

			const argv = [SCRIPT, '--bulk-csv', bulkCsv, '--out-dir', outDir, '--out-prefix', prefix,
				'--results-json', resultsJson];
			if (hasH5ad) {
				argv.push('--sc-h5ad', resolvePath(args.sc_h5ad, cwd));
			} else {
				argv.push('--sc-matrix-csv', resolvePath(args.sc_matrix_csv, cwd));
				argv.push('--sc-meta-csv', resolvePath(args.sc_meta_csv, cwd));
			}
			if (typeof args.group_col === 'string' && args.group_col.trim() !== '') {
				argv.push('--group-col', args.group_col.trim());
			}
			if (typeof args.gene_id_key === 'string' && args.gene_id_key.trim() !== '') {
				argv.push('--gene-id-key', args.gene_id_key.trim());
			}
			if (args.strip_version === true) argv.push('--strip-version');
			if (args.gene_id_case === 'upper') argv.push('--gene-to-upper');
			if (args.gene_id_case === 'lower') argv.push('--gene-to-lower');
			if (args.discriminability === true) argv.push('--discriminability');
			if (typeof args.normalize === 'string' && args.normalize !== '') {
				argv.push('--normalize', args.normalize);
			}
			if (typeof args.sc_lengths_col === 'string' && args.sc_lengths_col !== '') {
				argv.push('--sc-lengths-col', args.sc_lengths_col);
			}
			if (typeof args.bulk_lengths_col === 'string' && args.bulk_lengths_col !== '') {
				argv.push('--bulk-lengths-col', args.bulk_lengths_col);
			}
			if (args.no_manifest === true) argv.push('--no-manifest');
			if (args.no_figure === true) argv.push('--no-figure');

			const timeoutMs = Number.isFinite(args.timeout_ms) ? args.timeout_ms : defaultTimeoutMs;

			let outcome;
			try {
				outcome = await runScript(argv, { pythonPath, timeoutMs, logPath, cwd });
			} catch (error) {
				return {
					status: 'failed',
					message: `could not start the analysis interpreter (${pythonPath}): ${error.message}`,
					groups: [],
					warnings: [
						'Set pythonPath in this plugin row to an interpreter that has anndata, h5py, scipy, pandas, numpy and matplotlib.',
					],
					log_tail: [],
				};
			}

			const { lines, total } = await readTail(logPath, maxLogLines);

			if (outcome.signal !== null) {
				return {
					status: 'killed',
					message: `analysis killed after ${timeoutMs} ms (raise timeout_ms for large inputs)`,
					groups: [],
					warnings: [],
					log_tail: lines,
				};
			}

			let parsed;
			try {
				parsed = JSON.parse(await fs.readFile(resultsJson, 'utf8'));
			} catch {
				return {
					status: 'failed',
					message: `analysis exited with code ${outcome.code} and wrote no result file`,
					groups: [],
					warnings: [],
					log_tail: lines,
				};
			}

			if (parsed.schema !== RESULTS_SCHEMA) {
				return {
					status: 'failed',
					message: `unexpected result schema ${String(parsed.schema)}; expected ${RESULTS_SCHEMA}`,
					groups: [],
					warnings: [],
					log_tail: lines,
				};
			}

			const groups = Object.entries(parsed.pairs ?? {}).map(([group, s]) => ({
				group,
				n_cells: s.n_cells,
				n_genes_used: s.n_genes_used,
				pearson_r: s.pearson_r,
				spearman_rho: s.spearman_rho,
				slope: s.slope,
				intercept: s.intercept,
			}));

			const value = {
				status: 'ok',
				n_cells: parsed.scrna?.n_cells,
				n_genes_unique: parsed.scrna?.n_genes_unique,
				groups,
				warnings: [...(parsed.warnings ?? [])],
				results_json: resultsJson,
				manifest_json: parsed.outputs?.manifest_json ?? undefined,
				stats_csv: parsed.outputs?.stats_csv ?? undefined,
				figure_png: parsed.outputs?.figures?.[0],
				figure_pdf: parsed.outputs?.figures?.[1],
				log_tail: lines,
			};

			if (parsed.discriminability !== null && parsed.discriminability !== undefined) {
				const d = parsed.discriminability;
				value.discriminability = {
					matched_mean_r: d.matched_mean_r,
					mismatched_mean_r: d.mismatched_mean_r,
					gap: d.gap,
					n_groups: d.n_groups,
					n_correct: d.n_correct,
					misassigned: Object.entries(d.best_match ?? {})
						.filter(([, info]) => info.correct !== true)
						.map(([group, info]) => `${group}->${info.best_bulk_group}`),
				};
			}

			if (total > lines.length && lines.length > 0) {
				value.warnings.push(`log truncated to the last ${lines.length} of ${total} lines; full log at ${logPath}`);
			}

			return value;
		},
	};

	// Mounted on the plugin's own capability context, so unloading the plugin
	// removes the tool with it.
	ctx.effect(() => ctx.tools.register(tool));
}
