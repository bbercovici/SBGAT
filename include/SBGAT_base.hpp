#include <chrono>

class SBGAT_base{

public:
	SBGAT_base(){
		this -> modified();
	}
		

	void modified(){
		this -> time_modified = std::chrono::system_clock::now();
	}

protected:
	std::chrono::time_point<std::chrono::system_clock> time_modified;
	
};