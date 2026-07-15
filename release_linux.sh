#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>   (e.g. 0.1.0)"
    exit 1
fi

VERSION=$1
ROOT="$(cd "$(dirname "$0")" && pwd)"
RELEASE_DIR="$ROOT/releases"
NAME="ufbd-$VERSION-alpha-linux-x86_64"
STAGE="$RELEASE_DIR/$NAME"
BUILD="$ROOT/cmake-build-release"

rm -rf "$STAGE"
mkdir -p "$STAGE"

cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" > /dev/null
cmake --build "$BUILD" --target ufbd -j"$(nproc)"

cp "$BUILD/ufbd" "$STAGE/ufbd"
cp "$ROOT/README.md" "$STAGE/README.md"

cp "$BUILD/ufbd" "$ROOT/ufbd"

cd "$RELEASE_DIR"
rm -f "$NAME.zip"
zip -r "$NAME.zip" "$NAME" > /dev/null
rm -rf "$NAME"
echo "wrote $RELEASE_DIR/$NAME.zip ($(du -h "$NAME.zip" | cut -f1))"
