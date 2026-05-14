1) 'make' should compile the code
2) 'make clean' should clean executable and obj files

OBS: 
(i) In debug.c there is a function called "pause()"
if I compile in MAC, I need to use fpurge(stdin)
	and maybe comment out the line #include <stdio_ext.h> in the file header/utilities.h

if I compile in LINUX, I need to use __fpurge(stdin), 
	and maybe include the line #include <stdio_ext.h> in the file header/utilities.h

(ii) When compiling with LINUX, please load the mathematical library by including '-lm' in the LIBS variable, i.e.

#Libraries (Math)
LIBS= -lm


(if someone knows how to include if/else condition according to OS, it would be great....)