#!/bin/sh

type autoreconf &> /dev/null &&  { autoreconf --install; } || { 

{
	type libtoolize &>/dev/null &&
	{ echo 'running libtoolize'; libtoolize --copy --automake; } || 
	{ echo 'running glibtoolize'; glibtoolize --copy --automake; }
}  

{
	echo "running aclocal";
	aclocal;
} 

{
	echo "running autoheader [ignore the warnings]";
	autoheader;
}  

{
	echo "running automake";
	automake --foreign --add-missing --copy;
}  

{
	echo "running autoconf";
	autoconf;
}
} &&  
	echo "autogen complete"  ||
	echo "ERROR: autogen.sh failed, autogen is incomplete";

