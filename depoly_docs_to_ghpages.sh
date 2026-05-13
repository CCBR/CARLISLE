#!/usr/bin/env bash

set -euo pipefail

CONFIG_FILE="mkdocs.yml"

usage() {
  cat <<'EOF'
Usage:
  ./depoly_docs_to_ghpages.sh bootstrap <release-version>
  ./depoly_docs_to_ghpages.sh dev
  ./depoly_docs_to_ghpages.sh release <release-version>

Commands:
  bootstrap <release-version>
    First-time publication flow. Publishes:
      - dev
      - <release-version>
      - latest (alias to <release-version>)
    Then sets latest as the default docs version.

  dev
    Updates only the dev documentation.

  release <release-version>
    Publishes a new release version and moves latest to it.
    Does not change dev.

Examples:
  ./depoly_docs_to_ghpages.sh bootstrap 2.7.6
  ./depoly_docs_to_ghpages.sh dev
  ./depoly_docs_to_ghpages.sh release 2.7.7
EOF
}

require_version() {
  local version="${1:-}"
  if [[ -z "$version" ]]; then
    echo "Error: release version is required." >&2
    usage
    exit 1
  fi
}

deploy_dev() {
  mike deploy dev --push --config-file "$CONFIG_FILE"
}

deploy_release() {
  local version="$1"
  mike deploy "$version" latest --push --config-file "$CONFIG_FILE"
}

set_default_latest() {
  mike set-default latest --push
}

main() {
  local command="${1:-}"

  case "$command" in
    bootstrap)
      require_version "${2:-}"
      deploy_dev
      deploy_release "$2"
      set_default_latest
      ;;
    dev)
      deploy_dev
      ;;
    release)
      require_version "${2:-}"
      deploy_release "$2"
      set_default_latest
      ;;
    -h|--help|help|"")
      usage
      ;;
    *)
      echo "Error: unknown command '$command'." >&2
      usage
      exit 1
      ;;
  esac
}

main "$@"
