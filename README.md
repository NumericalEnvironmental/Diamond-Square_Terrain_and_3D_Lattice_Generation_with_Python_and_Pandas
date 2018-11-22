# Diamond-Square_Terrain_and_3D_Lattice_Generation_with_Python_and_Pandas

This is a relatively short python script (works with versions 2 or 3) for generating fractal surfaces or volumes, starting with an arbitrary initial lattice. By vectorizing operations using both numpy and pandas, the script can process fairly large lattices rapidly. Additional discussion and some example images are provided in my modeling and visualization blog (link pending).

The following two input files are required:

* lattice_input.csv - a comma-delimited text file with columns consisting of x-, y-, and z-coordinates, and an associated field value for each point. Lattices must be complete, with points sorted by x (fastest changing), then y, then z (last to change).
* params.txt - a text file listing the number of divisions along each dimension, in multiples of 2, along with a standard deviation to add random noise and a reduction factor for the standard deviation for each iteration.

The script will generate an output file (lattice_output.csv) formatted the same as the input file, with the resultant expanded lattice.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
