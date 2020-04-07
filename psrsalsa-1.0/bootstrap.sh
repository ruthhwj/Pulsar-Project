#/bin/sh

rm -rf autoconf
cp -r autoconf.boot autoconf

echo "Running autoconf..."

autoreconf --install --force

echo "If running autoconf was successful you should be able to run ./configure next."
