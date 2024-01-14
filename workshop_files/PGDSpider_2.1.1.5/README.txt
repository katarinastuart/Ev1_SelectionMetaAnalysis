..:: PGDSpider ::..

PGDSpider is a powerful automated data conversion tool for population genetic and genomics programs. 
It facilitates the data exchange possibilities between programs for a vast range of data types 
(e.g. DNA, RNA, NGS, microsatellite, SNP, RFLP, AFLP, multi-allelic data, allele frequency or genetic 
distances). Besides the conventional population genetics formats, PGDSpider integrates population 
genomics data formats commonly used to store and handle next-generation sequencing (NGS) data. 
Currently, PGDSpider is not meant to convert very large NGS files as it loads into memory the whole 
input file, whose size may exceed available RAM. However, since PGDSpider allows one to convert 
specific subsets of these NGS files into any other format, one could use this feature to calculate 
parameters or statistics for specific regions, and thus perform sliding window analysis over large 
genomic regions.

PGDSpider uses a newly developed PGD (Population Genetics Data) format as an intermediate step in the 
conversion process. PGD is a file format designed to store various kinds of population genetics data, 
including different data types (e.g. DNA sequences, microsatellites, AFLP or SNPs) and ploidy levels. 
PGD is based on the XML format and is therefore independent of any particular computer system and 
extensible for future needs. PGDSpider uses PGD to connect population genetics and genomics programs 
like a spider knits a web.

PGDSpider is written in Java and is therefore platform independent. It is user friendly due to its 
intuitive graphical user interface. PGDSpider allows the user to store his preferred conversion settings 
for repeated conversions of similar input formats. A command line version of PGDSpider is also provided, 
making it possible to embed PGDSpider in data analysis pipelines.



System requirements
====================
PGDSpider is written in Java and therefore platform independent, but SUN Java 1.7 RE (or a newer version) 
has to be installed. Java7 RE can be downloaded under following link: 
http://www.oracle.com/technetwork/java/javase/downloads/index.html



Installation Instructions
==========================

1st step:
----------

Install the Java7 RE

 * Windows: download and install Java7 RE with following link:
            http://www.oracle.com/technetwork/java/javase/downloads/index.html

 * Linux:   Ubuntu/Debian: Execute the following command as a root user:
                           "apt-get install openjdk-7-jre"

            Other LINUX distributions: http://www.oracle.com/technetwork/java/javase/downloads/index.html
 
 * Mac:     Apple Computer supplies their own version of Java. Use the Software Update feature (available on 
            the Apple menu) to check that you have the most up-to-date version of Java for your Mac. 
            Additionally, make sure that Java version 1.7 is set as first preference version. This can be 
            changed under "Applications - Utilities - Java Preferences.app". 
            If you have problems with downloading, installing or using Java on Mac, please contact Apple 
            Computer Technical Support.



2nd step:
----------
Unzip the PGDSpider_2.1.1.5.zip file on your local drive.

Execute PGDSpider GUI:

 * Windows:         execute the file "PGDSpider2.exe" to start the program
 
 * LINUX:           execute the command "./PGDSpider2.sh" to start the program   

 * Mac and others:  execute the command "java -Xmx1024m -Xms512m -jar PGDSpider2.jar" to start the program
 
 
Execute PGDSpider-cli (command line):

 * Windows:         execute the command "PGDSpider2-cli.exe"
 
 * Linux:           execute the command "java -Xmx1024m -Xms512M -jar PGDSpider2-cli.jar" 
 
 * Mac and others:  execute the command "java -Xmx1024m -Xms512M -jar PGDSpider2-cli.jar"



Java Web Start:
----------------
Additionally we provide the possibility to download and run PGDSpider from the web by the Java Web Start 
software. Java Web Start provides an easy, one-click activation of PGDSpider and it guarantees that you are 
always running the latest version.

Launch PGDSpider:
Java Web Start is included in the Java Runtime Environment. Have a look at the 1st step of the Installation 
Instructions to get information on how to get Java7 RE (or a newer version).

You can launch PGDSpider using Java Web Start from
 * Browser:           Click on the PGDSpider icon on the web page

 * Java Cache Viewer: To launch the PGDSpider Web Start a second time, you do not need to return to the web 
                      page where you first launched it; instead you can launch it from the Java Cache Viewer. 
                      To open the Java Cache Viewer execute following command in the console: javaws -viewer
                      To run PGDSpider Web Start, select it and click the Run button or double click the 
                      PGDSpider application.

 * Desktop:           You can add a desktop shortcut to the PGDSpider Web Start application. Select the 
                      application in the Java Cache Viewer (see above how to open it), then right-click and 
                      select “Install Shortcuts” or click the Install button. A shortcut is added to the 
                      desktop and you can launch the PGDSpider Web Start application just as you would launch 
                      any native application.

Limitations:
Starting PGDSpider from Java Web Start it is not possible to change the amount of memory PGDSpider is 
allowed to use (by default it is set to 1 GB). If you need to change the amount of memory (e.g.: if you
have large files to convert), download the PGDSpider application as described in the 2nd step of the 
Installation Instructions.



Help
=====
If you have any problems:
 * read the user manual delivered within the zip file
 * read the help file which can be found in the "Config" menu of the PGDSpider
 
 
 
How to cite PGDSpider and License
==================================
Lischer HEL and Excoffier L (2012) PGDSpider: An automated data conversion tool for connecting 
population genetics and genomics programs. Bioinformatics  28: 298-299.

Copyright (c) 2007-2018, Heidi E.L. Lischer. All rights reserved.

PGDSpider is distributed under the BSD 3-Clause License. 
For the full text of the license, see the file LICENSE.txt.
By using, modifying or distributing this software you agree to be bound by the terms of this license.



Contact and bug report
=======================
Heidi Lischer and Laurent Excoffier

Computational and Molecular Population Genetics lab (CMPG)
Institute of Ecology and Evolution (IEE)
University of Berne
3012 Bern
Switzerland

Members of the Swiss Institute of Bioinformatics (SIB)

e-mail: heidi.lischer@iee.unibe.ch


07. May 2018
