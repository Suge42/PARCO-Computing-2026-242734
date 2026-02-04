#!/bin/bash
set -e  # exit on first error

DESIRED_MATRIX=$1

if [ -z "$DESIRED_MATRIX" ]; then
  echo "Usage: $0 MATRIX_NAME"
  exit 1
fi

URL="https://suitesparse-collection-website.herokuapp.com/MM/HB/${DESIRED_MATRIX}.tar.gz"

wget -O "${DESIRED_MATRIX}.tar.gz" "$URL"
tar -xzf "${DESIRED_MATRIX}.tar.gz"

mv "${DESIRED_MATRIX}/${DESIRED_MATRIX}.mtx" "matrices/${DESIRED_MATRIX}.mtx"

rm -rf "${DESIRED_MATRIX}" "${DESIRED_MATRIX}.tar.gz"
