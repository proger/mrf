~This is the sample implementation of a Markov random field (MRF) based
image segmentation algorithm. The main code (mrf.cpp) has been written
by Csaba Gradwohl (Gradwohl.Csaba@stud.u-szeged.hu) with some minor
contributions from Zoltan Kato (kato@inf.u-szeged.hu,
http://www.inf.u-szeged.hu/~kato/). This code is released under the
GNU General Public License (see
http://www.gnu.org/copyleft/gpl.html). Please acknowledge the use of
our program by refering to the following papers (they contain detailed
information about the MRF model and optimization algorithms):

1) M. Berthod, Z. Kato, S. Yu, J. Zerubia: Bayesian image
classification using Markov random fields. Image and Vision Computing,
14(1996): 285-295, 1996.

2) Z. Kato: Multi-scale Markovian Modelisation in Computer Vision with
Applications to SPOT Image Segmentation. PhD thesis, INRIA Sophia
Antipolis, France, 1994.

3) Z. Kato, J. Zerubia and M. Berthod: Satellite image classification using a
modified Metropolis dynamics Proc. IEEE International Conf. on Acoust., Speech
and Sig. Proc., vol. 3, pp. 573-576, San Francisco, CA, March 23-26,
1992.

The program uses the "Mersenne Twister" random number generator
written by Agner Fog (http://www.agner.org/random/). The generator
itself is described in the article by M. Matsumoto & T. Nishimura, in:
ACM Transactions on Modeling and Computer Simulation, vol. 8, no. 1,
1998, pp. 3-30. Details on the initialization scheme can be found at
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html .

INSTALLATION:
=============

The code is platform independent. We have succesfully compiled it
under Linux (RedHat, Fedora, Debian) and Windows (XP, 2003 server).

1) You need wxWidgets 2.8 or later (http://www.wxwidgets.org)
2a) Under Linux/Unix, please edit the makefile under the "linux" subdirectory 
    and then
	$ make clean
	$ make
   should compile and install the program (it is called mrfdemo).
2b) Under Windows, we provide .NET 2005 (VC++8) compatible project
    files. Open the "mrfdemo" solution file under the "windows"
    subdirectory and choose "Build -> Rebuild mrfdemo" from the 
    menu to compile the program (it is called mrfdemo.exe). Make sure the 
    environment variable "WXWIN" is set correctly, otherwise you have to 
    modify the project settings. Note that only the "Release" configuration 
    works properly.
2b/1) Alternatively, you can also use the older wxWindows 2.4 compatible 
      solution files found under the subdirectory "windows-wxwin-2.4.x".
      We provide both VC++6 as well as .NET (VC++7) compatible project
      files. Of course, you need a properly installed wxWindows 2.4 in
      order to succesfully compile these versions.

USAGE NOTES:
============

The program works on BMP images. Some test images are provided under the
'images' subdirectory. The program GUI should be intuitive. Main
steps:

1) Load an image
2) Enter the number of pixel classes (~region type)
3) Push "Select classes" button
4) Press left mouse button over the input image and draw a rectangle
over a representative region of the first class. Then push "Next
class" button. The mean and variance should appear in the "Class
parameters" window. Continue with the next class until a
representative rectangle for all classes has been selected.
5) Set the weight of doubleton potentials (default is 0.9) and the
stopping threshold (iterations are stopped when the energy change is
less than the specified value).
6) Choose the optimization method from the pull-down list.
7) Adjust the optimization method's parameters:
	T0    - Initial temperature
	c     - temperature scheduler (T(n+1) = c*T(n))
	alpha - MMD's probability threshold
8) Push "Do it >>" button to execute segmentation.
9) Optionally, you can save the segmentation result as a BMP image.

During segmentation, the current classification along with the
temperature and global energy are displayed at each iteration. At the
end, the elapsed CPU time is also displayed (excluding GUI oveheaad!).

