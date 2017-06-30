# xyza2pipe

The xyza2pipe provides a cross conversion environment of higher dimensional NMR spectra (2D/3D/4D) between the following systems: NMRPipe, UCSF/Sparky, NMRView, XEASY (16bit)/CARA and Azara/ANSIG. It is possible to read VNMR (Agilent, formerly Varian) or XWinNMR (Bruker) binary spectra. This is open-source version branched from [**Olivia**](http://fermi.pharm.hokudai.ac.jp/olivia/) project.

## Conversion, Transposing, Combination and 2D Projection

- It supports transposing like NMRPipe does.
- Usage and parameter arguments are equivalent to NMRPipe.
- How to use:<br />
 Make a pipeline between importer and exporter of any combinations below.
- Conversion and transposing:<br />

Importer: File (reading) &rarr; Stream (stdout) | Exporter: Stream (stdin) &rarr; File (writing)
----------------------------------------------- | ----------------------------------------------
`xyza2pipe --in filename.ft` (2D)<br />`xyza2pipe --in filename%03d.ft` (3D)<br />`xyza2pipe --in filename%02d%03d.ft` (4D)<br />`ucsf2pipe --in filename.ucsf`<br />`nv2pipe --in filename.nv`<br />`xeasy2pipe --in filename.16`<br />`azara2pipe --in filename.spc`<br />`vnmr2pipe --in filename` (2D)<br />`vnmr2pipe --in filename%d` or `filename` (3D/4D)<br />`xwnmr2pipe --in 2rr/3rrr/4rrrr` | `pipe2xyza --out filename.ft` (2D)<br />`pipe2xyza --out filename%03d.ft` (3D)<br />`pipe2xyza --out filename%02d%03d.ft` (4D)<br />`pipe2ucsf --out filename.ucsf`<br />`pipe2nv --out filename.nv`<br />`pipe2xeasy --out filename.16`<br />`pipe2azara --out filename.spc`

- Combination:<br />

Importer: File (reading) &rarr; Stream (stdout) | Exporter: Stream (stdin) &rarr; File (writing)
----------------------------------------------- | ----------------------------------------------
`add2pipe --in1 file1.ft --in2 file2.ft`<br />`adducsf2pipe --in1 file1.ucsf --in2 file2.ucsf`<br />`addnv2pipe --in1 file1.nv --in2 file2.nv`<br />`addxeasy2pipe --in1 file1.16 --in2 file2.16`<br />`addazara2pipe --in1 file1.spc --in2 file2.spc`<br />`addvnmr2pipe --in1 file1 --in2 file2`<br />`addxwnmr2pipe --in1 file1 --in2 file2` | `pipe2xyza --out filename.ft` (2D)<br />`pipe2xyza --out filename%03d.ft` (3D)<br />`pipe2xyza --out filename%02d%03d.ft` (4D)<br />`pipe2ucsf --out filename.ucsf`<br />`pipe2nv --out filename.nv`<br />`pipe2xeasy --out filename.16`<br />`pipe2azara --out filename.spc`

- 2D Projection:<br />

Importer: File (reading) &rarr; Stream (stdout) | Exporter: Stream (stdin) &rarr; File (writing)
----------------------------------------------- | ----------------------------------------------
`xyza2pipe --in filename%03d.ft` (3D)<br />`xyza2pipe --in filename%02d%03d.ft` (4D)<br />`ucsf2pipe --in filename.ucsf`<br />`nv2pipe --in filename.nv`<br />`xeasy2pipe --in filename.16`<br />`azara2pipe --in filename.spc`<br />`vnmr2pipe --in filename%d` or `filename` (3D/4D)<br />`xwnmr2pipe --in 3rrr/4rrrr` | `pipe2proj --out proj.ft`<br />`pipe2proj | pipe2ucsf --out proj.ucsf`<br />`pipe2proj | pipe2nv --out proj.nv`<br />`pipe2proj | pipe2xeasy --out proj.16`<br />`pipe2proj | pipe2azara --out proj.spc`

- Effective pipeline combinations (conversion and transposing):<br />
 - To transpose NMRPipe formatted spectra, Type `xyza2pipe | pipe2xyza`.
 - To convert (UCSF&rarr;NMRPipe) and transpose, Type `ucsf2pipe | pipe2xyza`.
 - To convert (NMRView&rarr;NMRPipe) and transpose, Type `nv2pipe | pipe2xyza`.
 - To convert (XEASY&rarr;NMRPipe) and transpose, Type `xeasy2pipe | pipe2xyza`.
 - To convert (Azara&rarr;NMRPipe) and transpose, Type `azara2pipe | pipe2xyza`.<br /><br />
 - To transpose UCSF formatted spectra, Type `ucsf2pipe | pipe2ucsf`.
 - To convert (NMRView&rarr;UCSF) and transpose, Type `nv2pipe | pipe2ucsf`.
 - To convert (XEASY&rarr;UCSF) and transpose, Type `xeasy2pipe | pipe2ucsf`.
 - To convert (Azara&rarr;UCSF) and transpose, Type `azara2pipe | pipe2ucsf`.
 - To convert (NMRPipe&rarr;UCSF) and transpose, Type `xyza2pipe | pipe2ucsf`.<br /><br />
 - To transpose NMRView formatted spectra, Type `nv2pipe | pipe2nv`.
 - To convert (XEASY&rarr;NMRView) and transpose, Type `xeasy2pipe | pipe2nv`.
 - To convert (Azara&rarr;NMRView) and transpose, Type `azara2pipe | pipe2nv`.
 - To convert (NMRPipe&rarr;NMRView) and transpose, Type `xyza2pipe | pipe2nv`.
 - To convert (UCSF&rarr;NMRView) and transpose, Type `ucsf2pipe | pipe2nv`.<br /><br />
 - To transpose XEASY formatted spectra, Type `xeasy2pipe | pipe2xeasy`.
 - To convert (Azara&rarr;XEASY) and transpose, Type `azara2pipe | pipe2xeasy`.
 - To convert (NMRPipe&rarr;XEASY) and transpose, Type `xyza2pipe | pipe2xeasy`.
 - To convert (UCSF&rarr;XEASY) and transpose, Type `ucsf2pipe | pipe2xeasy`.
 - To convert (NMRView&rarr;XEASY) and transpose, Type `nv2pipe | pipe2xeasy`.<br /><br />
 - To transpose Azara formatted spectra, Type `azara2pipe | pipe2azara`.
 - To convert (NMRPipe&rarr;Azara) and transpose, Type `xyza2pipe | pipe2azara`.
 - To convert (UCSF&rarr;Azara) and transpose, Type `ucsf2pipe | pipe2azara`.
 - To convert (NMRView&rarr;Azara) and transpose, Type `nv2pipe | pipe2azara`.
 - To convert (XEASY&rarr;Azara) and transpose, Type `xeasy2pipe | pipe2azara`.<br /><br />
 - To convert (VNMR&rarr;NMRPipe) and transpose, Type `vnmr2pipe | pipe2xyza`.
 - To convert (VNMR&rarr;UCSF) and transpose, Type `vnmr2pipe | pipe2ucsf`.
 - To convert (VNMR&rarr;NMRView) and transpose, Type `vnmr2pipe | pipe2nv`.
 - To convert (VNMR&rarr;XEASY) and transpose, Type `vnmr2pipe | pipe2xeasy`.
 - To convert (VNMR&rarr;Azara) and transpose, Type `vnmr2pipe | pipe2azara`.<br /><br />
 - To convert (XWinNMR&rarr;NMRPipe) and transpose, Type `xwnmr2pipe | pipe2xyza`.
 - To convert (XWinNMR&rarr;UCSF) and transpose, Type `xwnmr2pipe | pipe2ucsf`.
 - To convert (XWinNMR&rarr;NMRView) and transpose, Type `xwnmr2pipe | pipe2nv`.
 - To convert (XWinNMR&rarr;XEASY) and transpose, Type `xwnmr2pipe | pipe2xeasy`.
 - To convert (XWinNMR&rarr;Azara) and transpose, Type `xwnmr2pipe | pipe2azara`.

- Example:<br />
 `ucsf2pipe --in ./hnco.ucsf | pipe2xyza --out ./data/hnco%03d.ft`<br />
 `nv2pipe   --in ./hnca.nv   | pipe2xyza --out ./data/hnca%03d.ft`

## Special notes

Monolithic NMRPipe format is compatible with this program. See examples below.
- Example (how to handle monolithic NMRPipe format):<br />
  `xyza2pipe --in  data%03d.ft > data.ft` [NMRPipe(regular)    &rarr; NMRPipe(monolithic)]<br />
  `pipe2xyza --out data%03d.ft < data.ft` [NMRPipe(monolithic) &rarr; NMRPipe(regular)]<br />
  `pipe2ucsf --out data.ucsf   < data.ft` [NMRPipe(monolithic) &rarr; UCSF]<br />
  `pipe2nv   --out data.nv     < data.ft` [NMRPipe(monolithic) &rarr; NMRView]

Several format conversions, such as NMRView&rarr;NMRPipe, UCSF&rarr;NMRPipe, VNMR&rarr;NMRPipe and XEASY&rarr;NMRPipe, cause to lose apodization information. This may affect result of NOESY back calculation in the **Olivia** system. For that cases, the conversion program assumes square of sine bell window function for direct observation axis (NMRPipe -fn SP -off 0.5 -end 0.98 -pow 2.0) and sine bell window function for the other axes (NMRPipe -fn SP -off 0.5 -end 0.98 -pow 1.0). If you want to retrieve correct apodization conditions, You should apply Fourier Transform again using either NMRPipe or Azara. On the other hand, Azara&rarr;NMRPipe conversion succeeds all apodization information and phase settings. While NMRPipe&rarr;Azara conversion doesn't.

**xeasy2pipe** requires a parameter file. The parameter file name should be *.param (parameter) for *.16 (XEASY spectra) where the asterisk (*) indicates common name.<br />
**azara2pipe** also requires a parameter file, which was implicitly generated by Azara at the execution time. The parameter file name should be *.spc.par and the asterisk (*) indicates common name of the Azara spectra *.spc.<br />
**vnmr2pipe** always requires a path to a directory that includes either procpar or procpar3d file.<br />
**xwnmr2pipe** requires properly situated parameter files, known as acqu*s and proc*s.

**pipe2proj** is used to make NMRPipe formatted 2D projection of spectra.

- Example:<br />
 `xyza2pipe --in ./data/hnco%03d.ft | pipe2proj --out hnco_proj.ft`

**add2pipe** is used to combine two NMRPipe formatted spectra like **addNMR** command.

- Example:<br />
 `add2pipe --in1 hsqc1.ft --in2 hsqc2.ft --add | pipe2xyza --out add_hsqc.ft`<br />
 `add2pipe --in1 hsqc1.ft --in2 hsqc2.ft --sub | pipe2xyza --out diff_hsqc.ft`

**defl2pipe** is used to deflate dimensionarity of pseudo arrayed NMRPipe 3D/4D spectra into regular 2D spectra.

## Resolving endianness

In relation to the [endianness](http://en.wikipedia.org/wiki/Endianness endianness) issue, program will try to detect proper endianess of file for CPU you used currently. However azara format is definitely endian dependent. If you fail in spectra conversion using **pipe2azara**, **pipe2nv**, **pipe2ucsf**, **pipe2xeasy** or **pipe2xyza**, Try to use '--swap' option which reverses byte order of data region.

**NOTE**: Improper endianness leads to **crash viewer applications**, to **generate unrecognized file format** or to **generate radial spectra**.

- Example (how to settle the endianness issue in the case of conversion from Azara to NMRPipe):<br />
 If the following pipeline command leads to generate radial spectra in some configuration.<br />
 `azara2pipe --in ./hsqc.spc | pipe2xyza --out ./hsqc.ft`
 
 In our experience, '--swap' option is helpful, where the option changes the byte order of data region of the NMRPipe file.<br />
 `azara2pipe --in ./hsqc.spc | pipe2xyza --swap --out ./hsqc.ft`

## Acknowledgments

**Bruce D. Ray**, NMR Center IUPUI, for general suggestions.<br />
**Miklos Guttman**, University of California, San Diego, for contribution to support UCSF format.<br />
**Alexander Eletsky**, University at Buffalo, for contribution to support XEASY format with code supply.<br />
**Tomohide Saio**, Rutgers University, for contribution to develop auto-repair mechanism against partially broken file.<br />
**Pascal Mercier**, University of Alberta, for alternative spectral regulation method relying on the origin frequency.<br />
**Naohiro Kobayashi**, Osaka University, for reporting NMRView spectra conversion issue.<br />
**Nicholas Fitzkee**, Mississippi State University, for reporting invalid Sparky spectra issue. 

## Release notes

**Oct 31, 2011**: ucsf2pipe and nv2pipe commands detect incorrect header information and fix them automatically.<br />
**Feb 14, 2012**: Fixed vnmr2pipe problems (distorted 2D spectra and shifted 3D spectra) and added support for varian 3D/4D complex spectra.<br />
**Feb 22, 2012**: Fixed Makefile.<br />
**Mar 16, 2012**: Added alternative spectral regulation method relying on the origin frequency. (--orig option)<br />
**Oct 01, 2012**: Fixed auto repairing routine to support pseudo-zerofill of ucsf and nv file.<br />
**Jan 30, 2014**: Added defl2pipe command to enforce dimensionarity to 2D spectra.<br />
**Jun 26, 2015**: Fixed invalid header of NMRView and UCSF spectra file.<br />
**Nov 16, 2015**: pipe2xyza allows inline transposing if filename is not specified.<br />
**Feb 07, 2017**: Released v1.0.0 @ GitHub, Inc. under Apache License V2.<br />
**Jun 30, 2017**: Released v1.0.1, fixed invalid header of UCSF spectra file.

## Executable files

### add2pipe: Two NMRPipe spectra &rarr; Combination &rarr; stdout stream
~~~
 add2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### addazara2pipe: Two Azara spectra &rarr; Combination &rarr; stdout stream
~~~
 addazara2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### addnv2pipe: Two NMRView spectra &rarr; Combination &rarr; stdout stream
~~~
 addnv2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### adducsf2pipe: Two UCSF spectra &rarr; Combination &rarr; stdout stream
~~~
 adducsf2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### adducsf2pipe: Two VNMR spectra &rarr; Combination &rarr; stdout stream
~~~
 addvnmr2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -p/--pdir1   A directory that includes VNMR procpar file1. (default=".")
  -q/--pdir2   A directory that includes VNMR procpar file2. (default=".")
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### addxeasy2pipe: Two XEASY spectra &rarr; Combination &rarr; stdout stream
~~~
 addxeasy2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### addxwnmr2pipe: Two XWinNMR spectra &rarr; Combination &rarr; stdout stream
~~~
 addxwnmr2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D)
  -a/--add inTemplate1 + inTemplate2 (default)
  -s/--sub inTemplate1 - inTemplate2
  -m/--mul inTemplate1 * inTemplate2
  --c1   Scale factor for inTemplate1 (default=1.0)
  --c2   Scale factor for inTemplate2 (default=1.0)
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### azara2pipe: Azara spectra &rarr; stdout stream
~~~
 azara2pipe -i/--in inTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### defl2pipe: Deflates dimensionarity into 2D spectra &rarr; stdout stream
~~~
 defl2pipe -i/--in inTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### nv2pipe: NMRView spectra &rarr; stdout stream
~~~
 nv2pipe -i/--in inTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
  --ndim  dimCount
~~~

### pipe2azara: stdin stream &rarr; Azara spectra
~~~
 pipe2azara -o/--out outTemplate(2D/3D/4D)
  -s/--swap  Byteswap data
  -l/--left  Left shift observable carrier frequency by sw/2
  -g/--orig  Rely on the origin frequency for spectral regulation
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### pipe2nv: stdin stream &rarr; NMRView spectra
~~~
 pipe2nv -o/--out outTemplate(2D/3D/4D)
  -s/--swap  Byteswap data
  -p/--pswap Byteswap parameter
  -l/--left  Left shift observable carrier frequency by sw/2
  -g/--orig  Rely on the origin frequency for spectral regulation
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### pipe2proj: stdin stream &rarr; Integration &rarr; NMRPipe spectra (with --out option) or stdout stream (without --out option)
~~~
 pipe2proj -o/--out outTemplate
  -s/--swap  Byteswap data
  -r/--right Right shift observable carrier frequency by sw/2
  -a/--abs   Absolute mode
  --xLAB X-Axis Label
  --yLAB Y-Axis Label
  --zLAB Z-Axis Label
  --aLAB A-Axis Label
  --xCAR X-Axis Center [ppm]
  --yCAR Y-Axis Center [ppm]
  --zCAR Z-Axis Center [ppm]
  --aCAR A-Axis Center [ppm]
~~~

### pipe2ucsf: stdin stream &rarr; UCSF spectra
~~~
 pipe2ucsf -o/--out outTemplate(2D/3D/4D)
  -s/--swap  Byteswap data
  -p/--pswap Byteswap parameter
  -l/--left  Left shift observable carrier frequency by sw/2
  -g/--orig  Rely on the origin frequency for spectral regulation
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### pipe2xeasy: stdin stream &rarr; XEASY spectra
~~~
 pipe2xeasy -o/--out outTemplate(2D/3D/4D)
  -s/--swap  Byteswap data
  -l/--left  Left shift observable carrier frequency by sw/2
  -g/--orig  Rely on the origin frequency for spectral regulation
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### pipe2xyza: stdin stream &rarr; NMRPipe spectra (with --out option) or stdout stream (without --out option)
~~~
 pipe2xyza -o/--out outTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  -s/--swap  Byteswap data
  -r/--right Right shift observable carrier frequency by sw/2
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### ucsf2pipe: UCSF spectra &rarr; stdout stream
~~~
 ucsf2pipe -i/--in inTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
  --ndim  dimCount
~~~

### vnmr2pipe: VNMR (Agilent, formerly Varian) binary spectra &rarr; stdout stream
~~~
 vnmr2pipe -i/--in inTemplate(2D/3D/4D)
  -p/--pdir    A directory that includes VNMR procpar file. (default=".")
  -e/--extLeft Extract left hand of X-Axis.
  -c/--adjCAR  Adjust indirect carrier frequencies according to 1H calibration.
  -h/--adjH2O  Adjust all carrier frequencies assuming transmitter is set to H2O solvent.
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
  --ndim  dimCount
  NOTE: User's phase correction is effective only in hyper-complex 2D spectra.
  --rp  0th Order X-Phase [deg]
  --lp  1st Order X-Phase [deg]
  --rp1 0th Order Y-Phase [deg]
  --lp1 1st Order Y-Phase [deg]
~~~

### xeasy2pipe: XEASY spectra &rarr; stdout stream
~~~
 xeasy2pipe -i/--in inTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~

### xwnmr2pipe: XWinNMR (Bruker) binary spectra &rarr; stdout stream
~~~
 xwnmr2pipe -i/--in inTemplate(2D/3D/4D)
  -e/--extLeft Extract left hand of X-Axis.
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
  --ndim  dimCount
~~~

### xyza2pipe: NMRPipe spectra &rarr; stdout stream
~~~
 xyza2pipe -i/--in inTemplate(2D/3D/4D)
  -x  Output X-Vector first
  -y  Output Y-Vector first
  -z  Output Z-Vector first
  -a  Output A-Vector first
  --xLAB  X-Axis Label
  --yLAB  Y-Axis Label
  --zLAB  Z-Axis Label
  --aLAB  A-Axis Label
  --xCAR  X-Axis Center [ppm]
  --yCAR  Y-Axis Center [ppm]
  --zCAR  Z-Axis Center [ppm]
  --aCAR  A-Axis Center [ppm]
~~~
