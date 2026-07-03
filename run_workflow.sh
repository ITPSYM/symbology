#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

mode="${1:-smoke}"
case "$mode" in
  smoke|full) ;;
  *)
    echo "Usage: $0 [smoke|full]" >&2
    exit 1
    ;;
esac

make

BOOTSTRAP="./bootstrap"
if [[ -x "./bootstrap.exe" ]]; then
  BOOTSTRAP="./bootstrap.exe"
fi

CONDITION="data/dlogmat_E6.wxf"
mkdir -p output logs

run_step() {
  local name="$1"
  shift
  echo "==> $name"
  "$@" 2>&1 | tee "logs/${name}.log"
}

extend_fec_to_6() {
  run_step extend_FEC_2 "$BOOTSTRAP" --extend -c "$CONDITION" -f data/FEC_1.wxf -o output/FEC_2.wxf
  run_step extend_FEC_3 "$BOOTSTRAP" --extend -c "$CONDITION" -f output/FEC_2.wxf -o output/FEC_3.wxf
  run_step extend_FEC_4 "$BOOTSTRAP" --extend -c "$CONDITION" -f output/FEC_3.wxf -o output/FEC_4.wxf
  run_step extend_FEC_5 "$BOOTSTRAP" --extend -c "$CONDITION" -f output/FEC_4.wxf -o output/FEC_5.wxf
  run_step extend_FEC_6 "$BOOTSTRAP" --extend -c "$CONDITION" -f output/FEC_5.wxf -o output/FEC_6.wxf
}

extend_lec_to_4() {
  run_step extend_LEC_2 "$BOOTSTRAP" --extend -c "$CONDITION" -l data/LEC_1.wxf -o output/LEC_2.wxf
  run_step extend_LEC_3 "$BOOTSTRAP" --extend -c "$CONDITION" -l output/LEC_2.wxf -o output/LEC_3.wxf
  run_step extend_LEC_4 "$BOOTSTRAP" --extend -c "$CONDITION" -l output/LEC_3.wxf -o output/LEC_4.wxf
}

if [[ "$mode" == "smoke" ]]; then
  extend_fec_to_6
  extend_lec_to_4
  run_step sew_2p2 "$BOOTSTRAP" --sew -c "$CONDITION" -f output/FEC_2.wxf -l output/LEC_2.wxf -o output/SEW_2p2.wxf
  run_step sew_3p1 "$BOOTSTRAP" --sew -c "$CONDITION" -f output/FEC_3.wxf -l data/LEC_1.wxf -o output/SEW_3p1.wxf
  run_step sew_4p2 "$BOOTSTRAP" --sew -c "$CONDITION" -f output/FEC_4.wxf -l output/LEC_2.wxf -o output/SEW_4p2.wxf
  run_step sew_5p1 "$BOOTSTRAP" --sew -c "$CONDITION" -f output/FEC_5.wxf -l data/LEC_1.wxf -o output/SEW_5p1.wxf
else
  extend_fec_to_6
  run_step extend_FEC_7 "$BOOTSTRAP" --extend -c "$CONDITION" -f output/FEC_6.wxf -o output/FEC_7.wxf
  run_step extend_FEC_8 "$BOOTSTRAP" --extend -c "$CONDITION" -f output/FEC_7.wxf -o output/FEC_8.wxf
  extend_lec_to_4
  run_step extend_LEC_5 "$BOOTSTRAP" --extend -c "$CONDITION" -l output/LEC_4.wxf -o output/LEC_5.wxf
  run_step sew_8p2 "$BOOTSTRAP" --sew -c "$CONDITION" -f output/FEC_8.wxf -l output/LEC_2.wxf -o output/SEW_8p2.wxf
fi
