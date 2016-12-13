#include <cmath>
#include "Tools.h"
#include "ImplicitQRSVD.h"
#include "SymmetricTridiagonal.h"
#include "SimulationDriver.h"
#include "EnergyTests.h"

void EnergyTest(){
  typedef double T;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
  int N=5;
  T a=(T)0,b=(T)1;
  T dX=(b-a)/(T)(N-1);
  JIXIE::NeoHookean<T> nh((T)1);
  JIXIE::LinearElasticity<T> le((T)1);
  JIXIE::FEMHyperelasticity<T> fem(a,dX,N,nh);
  TVect x(N);
  for(int i=0;i<N;i++) x(i)=(T).7*(a+dX*(T)i);

  JIXIE::EnergyTest<T> et("output",fem,10);
  et.RefinementTest(x);

}

void ElasticitySimulation(){


  typedef double T;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;

  JIXIE::ElasticityParameters<T> parameters;
  parameters.N=20;
  parameters.a=(T)0;
  T b=(T)1;

  parameters.c = (T)0.1; // added non-zero BC term

  parameters.dX=(b-parameters.a)/(T)(parameters.N-1);
  parameters.dt=(T).01;
  parameters.output_dir=std::string("output");
  parameters.rho=(T)1;
  parameters.k=(T)1;
  parameters.Newton_tol=(T)1e-8;
  parameters.max_newton_it=40;
  parameters.final_time=(T)4;
  parameters.frames_per_second=120;
  JIXIE::ElasticityDriver<T> driver(parameters);
  bool verbose=true;
  driver.RunSimulation(verbose);
}

void ConvertBinaryToDat(){
  typedef double T;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
  TVect x,v; int N=0,frame=0;

  std::string data_dir("output");
  std::string output_dat_dir("output/matlab");

  while(JIXIE::ElasticityDriver<T>::Read_State(x,v,N,data_dir,frame)){
    char str[12];
    sprintf(str, "%d", frame++);
    std::string frame_name(str);
    std::string positions_string(std::string("particle_x_")+frame_name);
    std::string velocities_string(std::string("particle_v_")+frame_name);
    FILE_IO::Write_DAT_File(std::string(output_dat_dir+std::string("/")+positions_string+std::string(".dat")),x);
    FILE_IO::Write_DAT_File(std::string(output_dat_dir+std::string("/")+velocities_string+std::string(".dat")),x);
  }
}

void ConvertBinaryToObj() {
	typedef double T;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVect;
	TVect x, v; int N = 0, frame = 0;

	std::string data_dir("output");
	std::string output_dat_dir("output/houdini");

	while (JIXIE::ElasticityDriver<T>::Read_State(x, v, N, data_dir, frame)) {
		char str[12];
		sprintf(str, "frame_%d", frame++);
		std::string frame_name(str);
		/*std::string positions_string(std::string("particle_x_") + frame_name);
		std::string velocities_string(std::string("particle_v_") + frame_name);
		FILE_IO::Write_DAT_File(std::string(output_dat_dir + std::string("/") + positions_string + std::string(".dat")), x);
		FILE_IO::Write_DAT_File(std::string(output_dat_dir + std::string("/") + velocities_string + std::string(".dat")), x);*/

		std::string file_name = std::string(output_dat_dir + std::string("/") + frame_name + std::string(".dat"));
		FILE* fpointer;
		fpointer = fopen(file_name.c_str(), "w");

		// print the vertices
		T dX = (T)1 / 19; // temporarily hardcoded
		T a0 = (T)1 /sqrt((x(1) - x(0)) / dX);
		// print first 4 nodes
		fprintf(fpointer, "v %f %f %f\n", x(0), a0, a0);
		fprintf(fpointer, "v %f %f %f\n", x(0), a0, -a0);
		fprintf(fpointer, "v %f %f %f\n", x(0), -a0, -a0);
		fprintf(fpointer, "v %f %f %f\n", x(0), -a0, a0);
		for (int i = 1; i < x.size() - 1; i++) {
			// compute a(x)
			T a = (T)1/sqrt((x(i + 1) - x(i - 1)) / 2 / dX);
			fprintf(fpointer, "v %f %f %f\n", x(i), a, a);
			fprintf(fpointer, "v %f %f %f\n", x(i), a, -a);
			fprintf(fpointer, "v %f %f %f\n", x(i), -a, -a);
			fprintf(fpointer, "v %f %f %f\n", x(i), -a, a);
		}
		// print last 4 nodes
		T aN= (T)1/sqrt((x(x.size()-1) - x(x.size()-2)) / dX);
		fprintf(fpointer, "v %f %f %f\n", x(x.size() - 1), aN, aN);
		fprintf(fpointer, "v %f %f %f\n", x(x.size() - 1), aN, -aN);
		fprintf(fpointer, "v %f %f %f\n", x(x.size() - 1), -aN, -aN);
		fprintf(fpointer, "v %f %f %f\n", x(x.size() - 1), -aN, aN);
		// print the triangles
		// first 2
		fprintf(fpointer, "f %i %i %i\n", 1, 2, 3);
		fprintf(fpointer, "f %i %i %i\n", 3, 4, 1);
		for (int i = 0; i < x.size(); i++) {
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 1, 4 * i + 2, 4 * i + 6);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 6, 4 * i + 5, 4 * i + 1);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 2, 4 * i + 3, 4 * i + 7);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 7, 4 * i + 6, 4 * i + 2);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 3, 4 * i + 4, 4 * i + 8);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 8, 4 * i + 7, 4 * i + 3);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 4, 4 * i + 1, 4 * i + 5);
			fprintf(fpointer, "f %i %i %i\n", 4 * i + 5, 4 * i + 8, 4 * i + 4);
		}// last 2
		fprintf(fpointer, "f %i %i %i\n", 4 * x.size() + 1, 4 * x.size() + 2, 4 * x.size() + 3);
		fprintf(fpointer, "f %i %i %i\n", 4 * x.size() + 3, 4 * x.size() + 4, 4 * x.size() + 1);

		fclose(fpointer);
	}
}

int main()
{
  //EnergyTest();
  //ElasticitySimulation();
  //ConvertBinaryToDat();
  ConvertBinaryToObj();
}
