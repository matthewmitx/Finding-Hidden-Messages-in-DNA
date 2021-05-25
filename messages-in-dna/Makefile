# A simple makefile for compiling three java classes
#

# define a makefile variable for the java compiler
#
JCC = javac
JVM = java

# define a makefile variable for compilation flags
# the -g flag compiles with debugging information
#
JFLAGS = -g
JFLAGSDIR = -d

MSGINDNA = .

# typing 'make' will invoke the first target entry in the makefile 
# (the default one in this case)
#
default: Helper.class ReverseComplement.class PatternMatching.class FrequentWords.class FindClumps.class MinimumSkew.class MotifEnumeration.class

main: Main.class

# this target entry builds the Average class
# the Average.class file is dependent on the Average.java file
# and the rule associated with this entry gives the command to create it

Helper.class: Helper.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) Helper.java

ReverseComplement.class: ReverseComplement.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) ReverseComplement.java

PatternMatching.class: PatternMatching.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) PatternMatching.java

FrequentWords.class: FrequentWords.java PatternMatching.java ReverseComplement.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) FrequentWords.java

FindClumps.class: FindClumps.java FrequentWords.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) FindClumps.java

MinimumSkew.class: MinimumSkew.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) MinimumSkew.java

MotifEnumeration.class: MotifEnumeration.java FrequentWords.java PatternMatching.java
	$(JCC) $(JFLAGSDIR) $(MSGINDNA) MotifEnumeration.java

Main.class: Main.java ReverseComplement.java PatternMatching.java FrequentWords.java FindClumps.java MinimumSkew.java MotifEnumeration.java
	$(JCC) $(JFLAGS) Main.java

# To start over from scratch, type 'make clean'.  
# Removes all .class files, so that the next make rebuilds them
#

clean:
	$(RM) *.class
	$(RM) MessagesInDNA/*.class