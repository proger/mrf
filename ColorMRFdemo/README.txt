This is the sample implementation of a Markov random field (MRF)
based color image segmentation algorithm. The main code (colormrf.cpp)
has been written by Mihaly Gara (gara@inf.u-szeged.hu,
http://www.inf.u-szeged.hu/~gara/) with some minor contributions from
Zoltan Kato (kato@inf.u-szeged.hu, http://www.inf.u-szeged.hu/~kato/)
using the intenisty-based segmentation code of Csaba Gradwohl. This
code is released under the GNU General Public License (see
http://www.gnu.org/copyleft/gpl.html). Please acknowledge the use of
our program by refering to the following paper:

1) Zoltan Kato, Ting Chuen Pong, and John Chung Mong Lee. Color Image
   Segmentation and Parameter Estimation in a Markovian
   Framework. Pattern Recognition Letters, 22(3-4):309--321, March 2001.

Note that the current demo program implements only a supervised
version of the segmentation method described in the above paper
(i.e. parameter values are learned interactively from representative
regions selected by the user). Otherwise, the program implements
exactly the color MRF model proposed in the paper.

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
   should compile and install the program (it is called colormrfdemo).
2b) Under Windows, we provide .NET 2005 (VC++8) compatible project
    files. Open the "colormrfdemo" solution file under the "windows"
    subdirectory and choose "Build -> Rebuild ColorMRFdemo" from the 
    menu to compile the program (it is called colormrfdemo.exe). Make sure the 
    environment variable "WXWIN" is set correctly, otherwise you have to 
    modify the project settings. Note that only the "Release" configuration 
    works properly.

USAGE NOTES:
============

The program works on BMP images. Some test images are provided under the
'images' subdirectory. The program GUI should be intuitive. Main
steps:

1) Load a color image
2) Enter the number of pixel classes (~region type)
3) Push "Select classes" button
4) Press left mouse button over the input image and draw a rectangle
over a representative region of the first class. Then push "Next
class" button. The mean and variance-covariance should appear in the "Class
parameters" window. Continue with the next class until a
representative rectangle for all classes has been selected.
5) Set the weight of doubleton potentials (default is 2.5) and the
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

