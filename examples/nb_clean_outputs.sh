#!/bin/bash

for file in "$@"; do
  tmp=$(mktemp)
  chmod --reference "$file" "$tmp"
  jq --indent 1 '
    ( .cells[] | select(.cell_type == "code") )
      |=
    ( .outputs = [] | .execution_count = null | del(.metadata.execution) )
    ' "$file" > "$tmp"
  mv "$tmp" "$file"
  # basename="${file%.*}"
  # extension="${file##*.}"
  # mv "$tmp" "${basename}.clean.${extension}"
done
