#!/usr/bin/env sh
# Print the body of one version's section from CHANGELOG.md.
# Usage: scripts/changelog_section.sh 0.1.0 [CHANGELOG.md]
# Only the first "## [VERSION]" heading is used; the section ends at the next
# "## " or "# " heading. Upstream OpenTFRaw history (after the first "# "
# heading below the title) is never searched. Exits non-zero if the section is
# missing or empty.
set -eu
version="$1"
file="${2:-CHANGELOG.md}"
body=$(awk -v heading="## [$version]" '
  found && (/^## / || /^# /) { exit }
  NR > 1 && /^# / { exit }
  found { print }
  !found && index($0, heading) == 1 { found = 1 }
' "$file")
if [ -z "$(printf '%s' "$body" | tr -d '[:space:]')" ]; then
  echo "error: no non-empty '## [$version]' section in $file" >&2
  exit 1
fi
printf '%s\n' "$body"
