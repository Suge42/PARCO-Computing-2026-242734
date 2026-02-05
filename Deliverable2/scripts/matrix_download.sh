#!/bin/bash
set -e  # exit on first error

MATRIX_GROUP=$1
DESIRED_MATRIX=$2

if [ -z "$DESIRED_MATRIX" ]; then
  echo "Usage: $0 MATRIX_NAME"
  exit 1
fi

URL="https://suitesparse-collection-website.herokuapp.com/MM/$MATRIX_GROUP/$DESIRED_MATRIX.tar.gz"

wget -O "$DESIRED_MATRIX.tar.gz" "$URL"
tar -xzf "$DESIRED_MATRIX.tar.gz"

mv "$DESIRED_MATRIX/$DESIRED_MATRIX.mtx" "matrices/$DESIRED_MATRIX.mtx"

rm -rf "$DESIRED_MATRIX" "$DESIRED_MATRIX.tar.gz"
