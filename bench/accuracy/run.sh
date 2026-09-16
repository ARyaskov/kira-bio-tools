#!/usr/bin/env bash
# Accuracy benchmark: call variants with the default pipeline and with each
# mpileup extension, then score every call set against a truth set.
#
#   bench/accuracy/run.sh <ref.fa> <reads.bam> <truth.vcf[.gz]> <confident.bed|-> <region|-> <outdir>
#
# Environment:
#   KIRA_BT   kira-bt binary (default target/release/kira-bt)
#   CONFIGS   space-separated "name=flags" list (default: every extension)
#   HAPPY     hap.py executable; when set, hap.py also scores each call set
#   MIN_QUAL  QUAL threshold applied to the calls (default 20)
#
# Writes <outdir>/<name>.vcf, <outdir>/<name>.json, <outdir>/results.tsv
# and <outdir>/results.md (a Markdown table for the CI job summary).
set -euo pipefail

ref=$1; bam=$2; truth=$3; bed=$4; region=$5; out=$6
KIRA_BT=${KIRA_BT:-target/release/kira-bt}
MIN_QUAL=${MIN_QUAL:-20}
here=$(cd "$(dirname "$0")" && pwd)
CONFIGS=${CONFIGS:-"default= hp-indel=--hp-indel indel-realign=--indel-realign_all assemble=--assemble recal=--recal nm-weight=--nm-weight_auto all=--hp-indel_--indel-realign_all_--assemble_--recal_--nm-weight_auto"}

mkdir -p "$out"
: > "$out/results.tsv"
region_arg=()
[ "$region" != "-" ] && region_arg=(-r "$region")
bed_arg=()
[ "$bed" != "-" ] && bed_arg=(--bed "$bed")
reg_eval=()
[ "$region" != "-" ] && reg_eval=(--region "$region")

for cfg in $CONFIGS; do
  name=${cfg%%=*}
  flags=${cfg#*=}
  flags=${flags//_/ }
  echo "== $name: mpileup $flags" >&2
  start=$(date +%s)
  # shellcheck disable=SC2086
  "$KIRA_BT" mpileup "$bam" $flags -o "$out/$name.mpileup.vcf" -- -f "$ref" "${region_arg[@]}"
  "$KIRA_BT" call "$out/$name.mpileup.vcf" -o "$out/$name.vcf" -- -mv
  secs=$(( $(date +%s) - start ))
  python3 "$here/evaluate.py" --truth "$truth" --query "$out/$name.vcf" "${bed_arg[@]}" "${reg_eval[@]}" \
    --min-qual "$MIN_QUAL" --label "$name" --json "$out/$name.json" --no-header \
    | sed "s/\$/\t$secs/" >> "$out/results.tsv"
  if [ -n "${HAPPY:-}" ]; then
    happy_args=()
    [ "$bed" != "-" ] && happy_args=(-f "$bed")
    [ "$region" != "-" ] && happy_args+=(-l "$region")
    "$HAPPY" "$truth" "$out/$name.vcf" -r "$ref" -o "$out/$name.happy" "${happy_args[@]}" >/dev/null 2>&1 \
      && cp "$out/$name.happy.summary.csv" "$out/$name.happy.csv" || echo "hap.py failed for $name" >&2
  fi
done

{
  echo "| config | class | truth | called | TP | FP | FN | precision | recall | F1 | GT conc. | seconds |"
  echo "|---|---|---|---|---|---|---|---|---|---|---|---|"
  awk -F'\t' '{ printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }' "$out/results.tsv"
} > "$out/results.md"
cat "$out/results.md"
