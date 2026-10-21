import argparse
import json
import re
from pathlib import Path

# This limit applies only to unstructured text previews.
# Structured tables are loaded in full.
MAX_PREVIEW_BYTES = 500_000


def data_lines(text):
    for line_number, line in enumerate(text.splitlines(), start=1):
        if line.strip() and not line.lstrip().startswith("#"):
            yield line_number, line


def table_result(item, headers, rows, errors=None):
    if errors:
        item["notice"] = (
            item.get("notice", "")
            + f" Skipped {len(errors)} malformed records. "
            + "Examples: "
            + "; ".join(errors[:3])
        ).strip()

    return {
        **item,
        "kind": "table",
        "headers": headers,
        "rows": rows,
    }


def parse_fixed_table(item, text, headers):
    rows = []
    errors = []

    for line_number, line in data_lines(text):
        fields = line.split()
        if len(fields) != len(headers):
            errors.append(f"line {line_number}: {line}")
        else:
            rows.append(fields)

    return table_result(item, headers, rows, errors)


def parse_named_text_table(item, text):
    lines = list(data_lines(text))

    if not lines:
        return table_result(item, ["Sequence", "MappedTo"], [])

    headers = lines[0][1].split()
    rows = []
    errors = []

    for line_number, line in lines[1:]:
        fields = line.split(maxsplit=len(headers) - 1)
        if len(fields) != len(headers):
            errors.append(f"line {line_number}: {line}")
        else:
            rows.append(fields)

    return table_result(item, headers, rows, errors)


def parse_check_table(item, text):
    headers = [
        "Sequence",
        "Query_Start",
        "Query_End",
        "Query_Coverage_Percent",
        "Mapped_to",
        "Mapped_Start",
        "Mapped_End",
        "Mapped_Start_Percent",
        "Mapped_End_Percent",
    ]

    pattern = re.compile(
        r"^(\S+)\s+(\d+)\s+(\d+)\s+"
        r"\(([\d.]+)\)\s+"
        r"(\S+)\s+(\d+)\s+(\d+)\s+"
        r"\(([\d.]+)\s*,\s*([\d.]+)\)\s*$"
    )

    rows = []
    errors = []

    for line_number, line in data_lines(text):
        match = pattern.match(line.strip())
        if match:
            rows.append(list(match.groups()))
        else:
            errors.append(f"line {line_number}: {line}")

    return table_result(item, headers, rows, errors)


def parse_target_attributes(attributes):
    # Split on semicolons outside quoted annotation text.
    parts = re.split(
        r';(?=(?:[^"]*"[^"]*")*[^"]*$)',
        attributes,
    )

    values = {}
    for part in parts:
        match = re.match(
            r"^\s*([A-Za-z][A-Za-z0-9_]*)"
            r"(?:\s*=\s*|\s+)(.*?)\s*$",
            part,
        )
        if match:
            values[match.group(1)] = match.group(2)

    target = values.get("Target", "").split()
    if len(target) < 3:
        raise ValueError("Missing or malformed Target attribute")

    return [
        target[0],
        target[1],
        target[2],
        values.get("QueryLength", values.get("TargetQueryLength", "")),
        values.get("TargetLength", ""),
        values.get("Annot", ""),
    ]


def parse_gff(item, text, stage):
    base_headers = [
        "Contig", "Source", "Type", "Start",
        "End", "Score", "Strand", "Phase",
    ]

    query_headers = [
        "Contig", "Source", "Query_Type", "Query_Start",
        "Query_End", "Query_Score", "Strand", "Query_Phase",
    ]

    filter_stages = {"07_filter_sequences", "07_filter_sequence"}

    if stage == "03_seqer":
        headers = base_headers + ["Family", "Length"]

    elif stage == "06_classify_consensus":
        headers = query_headers + [
            "Length", "Type", "Mapped_to",
            "Phase", "Start", "End", "Score",
        ]
        item["notice"] = (
            "Query_* columns contain the original GFF fields. "
            "The expanded Score column follows the requested label. "
            "In the current CARP generator, this attribute represents "
            "the remaining reference length, not an alignment score."
        )

    elif stage in filter_stages:
        headers = query_headers + [
            "Length", "Target", "Start", "End",
            "QueryLength", "TargetLength", "Annot",
        ]
        item["notice"] = (
            "Query_* columns contain the original GFF fields. "
            "Start and End contain Target coordinates."
        )

    else:
        # Do not add an artificial Family column.
        headers = base_headers + ["Attributes", "Length"]

    rows = []
    errors = []

    for line_number, line in enumerate(text.splitlines(), start=1):
        if line.strip() == "##FASTA":
            break
        if not line.strip() or line.lstrip().startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) != 9:
            errors.append(f"line {line_number}: expected nine GFF fields")
            continue

        try:
            length = int(fields[4]) - int(fields[3]) + 1
            attributes = fields[8]

            if stage == "06_classify_consensus":
                mapped = attributes.split()
                if len(mapped) != 6:
                    raise ValueError("Expected six mapping attribute fields")
                row = fields[:8] + [length] + mapped

            elif stage in filter_stages:
                row = (
                    fields[:8]
                    + [length]
                    + parse_target_attributes(attributes)
                )

            else:
                row = fields[:8] + [attributes, length]

            rows.append(row)

        except ValueError as error:
            errors.append(f"line {line_number}: {error}")

    return table_result(item, headers, rows, errors)


def parse_repeatmasker_out(item, text):
    headers = [
        "SW_Score",
        "Divergence_Percent",
        "Deletion_Percent",
        "Insertion_Percent",
        "Query",
        "Query_Start",
        "Query_End",
        "Query_Left",
        "Strand",
        "Repeat",
        "Class_Family",
        "Repeat_Start",
        "Repeat_End",
        "Repeat_Left",
        "ID",
        "Overlap",
        "Query_Length",
    ]

    rows = []
    errors = []

    def integer(value):
        return int(value.strip("()"))

    for line_number, line in data_lines(text):
        fields = line.split()

        # Skip standard headers and status messages.
        if not fields or not fields[0].isdigit():
            continue

        try:
            if len(fields) not in {15, 16}:
                raise ValueError("Expected 15 fields and an optional overlap flag")

            if len(fields) == 16 and fields[15] != "*":
                raise ValueError("Unrecognized overlap flag")

            strand = fields[8]

            if strand == "C":
                # Reverse-strand layout: (left), end, start.
                repeat_left = integer(fields[11])
                repeat_end = integer(fields[12])
                repeat_start = integer(fields[13])
            elif strand == "+":
                repeat_start = integer(fields[11])
                repeat_end = integer(fields[12])
                repeat_left = integer(fields[13])
            else:
                raise ValueError(f"Unrecognized strand: {strand}")

            query_start = integer(fields[5])
            query_end = integer(fields[6])

            rows.append([
                int(fields[0]),
                float(fields[1]),
                float(fields[2]),
                float(fields[3]),
                fields[4],
                query_start,
                query_end,
                integer(fields[7]),
                strand,
                fields[9],
                fields[10],
                repeat_start,
                repeat_end,
                repeat_left,
                integer(fields[14]),
                "*" if len(fields) == 16 else "",
                query_end - query_start + 1,
            ])

        except ValueError as error:
            errors.append(f"line {line_number}: {error}")

    return table_result(item, headers, rows, errors)


def parse_tsv(item, text):
    lines = list(data_lines(text))
    if not lines:
        return table_result(item, ["Value"], [])

    headers = lines[0][1].split("\t")
    rows = []
    errors = []

    for line_number, line in lines[1:]:
        fields = line.split("\t")
        if len(fields) == len(headers):
            rows.append(fields)
        else:
            errors.append(f"line {line_number}: unexpected column count")

    return table_result(item, headers, rows, errors)


def read_file(path, root):
    item = {
        "name": str(path.relative_to(root)),
        "notice": "",
    }

    stage = root.name
    suffix = path.suffix.lower()
    name = path.name.lower()

    is_gff = suffix in {".gff", ".gff3"}
    is_index = stage == "04_finalize_consensus" and suffix == ".idx"
    is_map = stage == "05_repeats_in_repbase" and suffix == ".map"
    is_classification = stage == "06_classify_consensus" and suffix == ".txt"
    is_protein = (
        stage == "08_get_protein_info"
        and name.endswith("_protein.txt")
    )
    is_repeatmasker = (
        stage == "12_genome_annotation"
        and name.endswith(".mfa.out")
    )

    structured = (
        is_gff or is_index or is_map or is_classification
        or is_protein or is_repeatmasker or suffix == ".tsv"
    )

    if structured:
        # Read complete tables so filtering uses every parsed record.
        text = path.read_text(encoding="utf-8", errors="replace")

        if is_gff:
            return parse_gff(item, text, stage)

        if is_index:
            item["notice"] = (
                "Index_Value_1 and Index_Value_2 preserve the original "
                "numeric fields without assigning unverified meanings."
            )
            return parse_fixed_table(
                item, text,
                ["File", "Sequence", "Index_Value_1", "Index_Value_2"],
            )

        if is_map:
            item["notice"] = (
                "Field_8, Field_9, Field_11, and Field_12 retain "
                "their original positions because their meanings "
                "have not been verified for this output."
            )
            return parse_fixed_table(
                item, text,
                [
                    "Sequence", "Query_Start", "Query_End",
                    "Mapped_to", "Mapped_Start", "Mapped_End",
                    "Orientation", "Field_8", "Field_9",
                    "Score", "Field_11", "Field_12",
                ],
            )

        if is_classification:
            if name.endswith("_check.txt"):
                return parse_check_table(item, text)
            return parse_named_text_table(item, text)

        if is_protein:
            return parse_named_text_table(item, text)

        if is_repeatmasker:
            return parse_repeatmasker_out(item, text)

        return parse_tsv(item, text)

    with path.open("rb") as handle:
        raw = handle.read(MAX_PREVIEW_BYTES + 1)

    truncated = len(raw) > MAX_PREVIEW_BYTES
    raw = raw[:MAX_PREVIEW_BYTES]

    if b"\x00" in raw:
        return {
            **item,
            "kind": "text",
            "text": "Binary file. Preview is unavailable.",
        }

    text = raw.decode("utf-8", errors="replace")
    if truncated:
        text = text.rsplit("\n", 1)[0] if "\n" in text else ""
        item["notice"] = (
            f"This file exceeds {MAX_PREVIEW_BYTES:,} bytes. "
            "Only a partial text preview is displayed."
        )

    return {**item, "kind": "text", "text": text}


PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CARP Output Report</title>
<style>
* { box-sizing: border-box; }
body {
  margin: 0;
  font: 14px system-ui, sans-serif;
  color: #243247;
  background: #f5f7fb;
}
.layout { display: flex; min-height: 100vh; }
aside {
  width: 280px;
  flex-shrink: 0;
  padding: 24px 14px;
  background: #162238;
  color: white;
}
aside h2 { margin: 0 8px 8px; }
aside p { margin: 0 8px 24px; color: #b9c6db; }
.tab {
  display: block;
  width: 100%;
  text-align: left;
  overflow-wrap: anywhere;
  margin: 7px 0;
  padding: 12px;
  border: 0;
  border-radius: 7px;
  background: transparent;
  color: #dbe5f6;
  cursor: pointer;
}
.tab:hover { background: #263751; }
.tab.active { background: #2563eb; color: white; }
main { flex: 1; min-width: 0; padding: 28px; }
h1 { margin-top: 0; overflow-wrap: anywhere; }
.card {
  background: white;
  border: 1px solid #e0e6ef;
  border-radius: 10px;
  padding: 18px;
}
.controls {
  display: flex;
  gap: 10px;
  flex-wrap: wrap;
  align-items: center;
  margin: 14px 0;
}
input, select, button {
  font: inherit;
  padding: 8px 10px;
  border: 1px solid #cdd6e3;
  border-radius: 5px;
}
button { background: white; cursor: pointer; }
button:hover { background: #eef3ff; }
button:disabled { opacity: .45; cursor: default; }
input[type=number] { width: 145px; }
#files { max-width: 100%; }
#notice {
  color: #975b00;
  white-space: pre-wrap;
  overflow-wrap: anywhere;
}
#status { color: #5c6c82; margin: 12px 0; }
.help { color: #5c6c82; line-height: 1.6; }
.table-wrap { overflow: auto; max-height: 68vh; }
table { width: 100%; border-collapse: collapse; white-space: nowrap; }
th, td {
  text-align: left;
  padding: 9px 12px;
  border-bottom: 1px solid #e8edf4;
}
th { position: sticky; top: 0; background: #eef3fb; }
th button {
  padding: 0;
  border: 0;
  background: transparent;
  font-weight: 600;
}
tr:nth-child(even) { background: #fafcff; }
pre {
  overflow: auto;
  max-height: 72vh;
  padding: 16px;
  background: #f7f9fc;
  white-space: pre;
}
[hidden] { display: none !important; }
@media (max-width: 760px) {
  .layout { flex-direction: column; }
  aside { width: 100%; }
  main { padding: 16px; }
}
</style>
</head>
<body>
<div class="layout">
  <aside>
    <h2>CARP-snakemake report</h2>
    <p id="sample"></p>
    <nav id="tabs" aria-label="Output folders"></nav>
  </aside>

  <main>
    <h1 id="folderTitle"></h1>
    <div class="card">
      <label>File: <select id="files"></select></label>
      <p id="notice"></p>

      <section id="tablePanel" hidden>
        <div class="controls">
          <input id="search" placeholder="Search all columns"
                 aria-label="Search all columns">
          <select id="column" aria-label="Numeric column"></select>
          <input id="lower" type="number" step="any"
                 placeholder="Minimum value" aria-label="Minimum value">
          <input id="upper" type="number" step="any"
                 placeholder="Maximum value" aria-label="Maximum value">
          <button id="minimum">Show minimum</button>
          <button id="maximum">Show maximum</button>
          <button id="reset">Reset filters</button>
        </div>

        <p class="help">
          Click a column header to sort.
          Range limits are inclusive.
          Minimum and maximum values are calculated after applying
          search and range filters. All tied records are retained.
          Missing numeric values are excluded from numeric filtering.
        </p>

        <div id="status" role="status"></div>

        <div class="table-wrap">
          <table>
            <thead id="thead"></thead>
            <tbody id="tbody"></tbody>
          </table>
        </div>

        <div class="controls">
          <button id="previous">Previous</button>
          <span id="pageInfo"></span>
          <button id="next">Next</button>
        </div>
      </section>

      <pre id="textPanel" hidden></pre>
    </div>
  </main>
</div>

<script id="reportData" type="application/json">__REPORT_DATA__</script>
<script>
"use strict";

const DATA = JSON.parse(document.getElementById("reportData").textContent);
const $ = id => document.getElementById(id);
const PAGE_SIZE = 100;

let folderIndex = 0;
let file = null;
let numericColumns = [];
let sortColumn = -1;
let ascending = true;
let page = 0;
let extreme = null;

function number(value) {
  if (value === null || value === undefined) return null;
  const text = String(value).trim();
  if (text === "" || text === ".") return null;
  const result = Number(text);
  return Number.isFinite(result) ? result : null;
}

function isNumericColumn(index) {
  let found = false;
  for (const row of file.rows) {
    const value = row[index];
    if (value == null || String(value).trim() === "" || value === ".") {
      continue;
    }
    if (number(value) === null) return false;
    found = true;
  }
  return found;
}

function option(text, value) {
  const element = document.createElement("option");
  element.textContent = text;
  element.value = value;
  return element;
}

function selectFolder(index) {
  folderIndex = index;

  document.querySelectorAll(".tab").forEach((button, i) => {
    button.classList.toggle("active", i === index);
    button.setAttribute("aria-pressed", String(i === index));
  });

  const folder = DATA.folders[index];
  $("folderTitle").textContent = folder.name;
  $("files").replaceChildren();

  folder.files.forEach((entry, i) => {
    $("files").append(option(entry.name, i));
  });

  $("files").disabled = folder.files.length === 0;

  if (folder.files.length) {
    selectFile(0);
  } else {
    $("notice").textContent = "No files are available in this folder.";
    $("tablePanel").hidden = true;
    $("textPanel").hidden = true;
  }
}

function clearFilters() {
  $("search").value = "";
  $("lower").value = "";
  $("upper").value = "";
  extreme = null;
  page = 0;
}

function selectFile(index) {
  file = DATA.folders[folderIndex].files[index];
  $("notice").textContent = file.notice || "";
  $("tablePanel").hidden = file.kind !== "table";
  $("textPanel").hidden = file.kind !== "text";

  if (file.kind === "text") {
    $("textPanel").textContent = file.text;
    return;
  }

  clearFilters();
  sortColumn = -1;
  ascending = true;

  numericColumns = file.headers.map((_, index) => isNumericColumn(index));
  $("column").replaceChildren(option("Select numeric column", ""));

  file.headers.forEach((header, index) => {
    if (numericColumns[index]) {
      $("column").append(option(header, index));
    }
  });

  const lengthIndex = file.headers.indexOf("Length");
  const firstNumeric = numericColumns.indexOf(true);
  const defaultColumn = (
    lengthIndex >= 0 && numericColumns[lengthIndex]
      ? lengthIndex
      : firstNumeric
  );

  if (defaultColumn >= 0) {
    $("column").value = String(defaultColumn);
  }

  render();
}

function render() {
  let rows = file.rows.slice();
  const query = $("search").value.trim().toLowerCase();
  const selected = $("column").value;
  const column = selected === "" ? null : Number(selected);
  const lower = number($("lower").value);
  const upper = number($("upper").value);

  $("minimum").disabled = column === null;
  $("maximum").disabled = column === null;
  $("lower").disabled = column === null;
  $("upper").disabled = column === null;

  if (query) {
    rows = rows.filter(row =>
      row.some(value => String(value).toLowerCase().includes(query))
    );
  }

  if (column !== null && (lower !== null || upper !== null)) {
    rows = rows.filter(row => {
      const value = number(row[column]);
      return value !== null
        && (lower === null || value >= lower)
        && (upper === null || value <= upper);
    });
  }

  let extremeText = "";

  if (extreme && column !== null) {
    let target = null;

    for (const row of rows) {
      const value = number(row[column]);
      if (value === null) continue;
      if (
        target === null ||
        (extreme === "min" && value < target) ||
        (extreme === "max" && value > target)
      ) {
        target = value;
      }
    }

    if (target === null) {
      rows = [];
      extremeText = " | No valid numeric values";
    } else {
      rows = rows.filter(row => number(row[column]) === target);
      extremeText = ` | ${file.headers[column]} ${
        extreme === "min" ? "minimum" : "maximum"
      } = ${target}`;
    }
  }

  if (sortColumn >= 0) {
    const numeric = numericColumns[sortColumn];

    rows.sort((a, b) => {
      if (numeric) {
        const x = number(a[sortColumn]);
        const y = number(b[sortColumn]);

        // Keep missing values last in either sorting direction.
        if (x === null) return y === null ? 0 : 1;
        if (y === null) return -1;
        return ascending ? x - y : y - x;
      }

      const order = String(a[sortColumn]).localeCompare(
        String(b[sortColumn]), "en", {numeric: true}
      );
      return ascending ? order : -order;
    });
  }

  const pages = Math.max(1, Math.ceil(rows.length / PAGE_SIZE));
  page = Math.max(0, Math.min(page, pages - 1));

  const invalidRange = (
    column !== null && lower !== null && upper !== null && lower > upper
  );

  $("status").textContent =
    `Loaded: ${file.rows.length} records | Matching: ${rows.length}` +
    extremeText +
    (invalidRange ? " | Minimum exceeds maximum" : "");

  $("pageInfo").textContent = `Page ${page + 1} of ${pages}`;
  $("previous").disabled = page === 0;
  $("next").disabled = page >= pages - 1;

  const headerRow = document.createElement("tr");

  file.headers.forEach((header, index) => {
    const th = document.createElement("th");
    const button = document.createElement("button");

    th.setAttribute(
      "aria-sort",
      sortColumn === index
        ? (ascending ? "ascending" : "descending")
        : "none"
    );

    button.textContent = header + (
      sortColumn === index ? (ascending ? " ▲" : " ▼") : ""
    );

    button.onclick = () => {
      ascending = sortColumn === index ? !ascending : true;
      sortColumn = index;
      page = 0;
      render();
    };

    th.append(button);
    headerRow.append(th);
  });

  $("thead").replaceChildren(headerRow);

  const fragment = document.createDocumentFragment();
  const visibleRows = rows.slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE);

  for (const row of visibleRows) {
    const tr = document.createElement("tr");
    for (const value of row) {
      const td = document.createElement("td");
      td.textContent = value;
      tr.append(td);
    }
    fragment.append(tr);
  }

  if (!visibleRows.length) {
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = file.headers.length;
    td.textContent = "No matching records.";
    tr.append(td);
    fragment.append(tr);
  }

  $("tbody").replaceChildren(fragment);
}

$("files").onchange = () => selectFile(Number($("files").value));

["search", "lower", "upper"].forEach(id => {
  $(id).addEventListener("input", () => {
    extreme = null;
    page = 0;
    render();
  });
});

$("column").onchange = () => {
  $("lower").value = "";
  $("upper").value = "";
  extreme = null;
  page = 0;
  render();
};

$("minimum").onclick = () => {
  extreme = "min";
  page = 0;
  render();
};

$("maximum").onclick = () => {
  extreme = "max";
  page = 0;
  render();
};

$("reset").onclick = () => {
  clearFilters();
  render();
};

$("previous").onclick = () => {
  page--;
  render();
};

$("next").onclick = () => {
  page++;
  render();
};

$("sample").textContent = DATA.sample;

DATA.folders.forEach((folder, index) => {
  const button = document.createElement("button");
  button.className = "tab";
  button.textContent = folder.name;
  button.onclick = () => selectFolder(index);
  $("tabs").append(button);
});

if (DATA.folders.length) {
  selectFolder(0);
} else {
  $("folderTitle").textContent = "No output folders found";
  $("files").disabled = true;
}
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(
        description="Create an interactive HTML report from CARP output folders."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Sample output directory, such as output/test_input",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("report.html"),
        help="Destination HTML file",
    )
    args = parser.parse_args()

    root = args.input_dir.resolve()
    destination = args.output.resolve()

    if not root.is_dir():
        parser.error(f"Directory does not exist: {root}")

    folders = []

    for folder in sorted(root.iterdir()):
        if (
            not folder.is_dir()
            or folder.is_symlink()
            or folder.name.startswith(".")
        ):
            continue

        files = []

        for path in sorted(folder.rglob("*")):
            if not path.is_file() or path.is_symlink():
                continue
            if path.resolve() == destination:
                continue
            if any(
                part.startswith(".")
                for part in path.relative_to(folder).parts
            ):
                continue

            try:
                entry = read_file(path, folder)
            except OSError as error:
                entry = {
                    "name": str(path.relative_to(folder)),
                    "kind": "text",
                    "notice": "",
                    "text": f"Unable to read file: {error}",
                }

            files.append(entry)

        folders.append({
            "name": folder.name,
            "files": files,
        })

    data = json.dumps(
        {"sample": root.name, "folders": folders},
        ensure_ascii=False,
        allow_nan=False,
    )

    # Prevent file content from closing the embedded JSON script element.
    data = (
        data.replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        PAGE.replace("__REPORT_DATA__", data),
        encoding="utf-8",
    )

    file_count = sum(len(folder["files"]) for folder in folders)
    print(f"Included {len(folders)} folders and {file_count} files.")
    print(f"Report generated: {destination}")


if __name__ == "__main__":
    main()