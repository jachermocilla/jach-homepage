#!/bin/bash
#jach

echo "Preparing...."
rm /etc/named.conf 2>&1 > /dev/null
rm /etc/named.run 2>&1 > /dev/null
rm -fr /var/named 2>&1 > /dev/null
mkdir /var/named 2>&1 > /dev/null
chown root.root /var/named 2>&1 > /dev/null
chown root.root /var/run/named 2>&1 > /dev/null
echo "done."

