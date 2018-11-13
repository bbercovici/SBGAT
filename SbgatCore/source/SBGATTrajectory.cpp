#include <SBGATTrajectory.hpp>
#include <OrbitConversions.hpp>



void SBGATTrajectory::GenerateKeplerianTrajectory(std::vector<arma::vec> & positions,
	std::vector<arma::vec> & velocities,
	const std::vector<double> & times,const arma::vec & elements,const double & mu){

	OC::KepState kep(elements,mu);
	OC::CartState state;

	for (const auto time : times){
		state = kep.convert_to_cart(time);

		positions.push_back(state.get_position_vector());
		velocities.push_back(state.get_velocity_vector());

	}

}
