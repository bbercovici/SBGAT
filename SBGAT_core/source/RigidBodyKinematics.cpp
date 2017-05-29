#include "RigidBodyKinematics.hpp"
#include <math.h>


arma::mat mrp_to_dcm(const arma::vec & sigma) {
	arma::mat identity(3, 3);
	identity = identity.eye();

	arma::mat dcm = identity + (8 * tilde(sigma) * (tilde(sigma))
	                            - 4 * (1 - pow(arma::norm(sigma), 2)) * tilde(sigma)) / pow(1 + pow(arma::norm(sigma), 2), 2);
	return dcm;
}


arma::mat euler321_to_dcm(const arma::vec & euler_angles) {
	arma::mat M(M1(euler_angles(2)) * M2(euler_angles(1)) * M3(euler_angles(0)));
	return M;
}


arma::mat longitude_latitude_to_dcm(const arma::vec & euler_angles) {
	arma::mat M(M2(-euler_angles(1)) * M3(euler_angles(0)));
	return M;
}


arma::mat euler313_to_dcm(const arma::vec & euler_angles) {
	arma::mat M(M3(euler_angles(2)) * M1(euler_angles(1)) * M3(euler_angles(0)));
	return M;
}

arma::mat euler321d_to_dcm(arma::vec & euler_angles) {
	euler_angles = arma::datum::pi / 180. * euler_angles;
	arma::mat M(M1(euler_angles(2)) * M2(euler_angles(1)) * M3(euler_angles(0)));
	return M;
}

arma::mat euler313d_to_dcm(arma::vec & euler_angles) {
	euler_angles = arma::datum::pi / 180. * euler_angles;

	arma::mat M(M3(euler_angles(2)) * M1(euler_angles(1)) * M3(euler_angles(0)));
	return M;
}

arma::vec mrp_to_quat(const arma::vec & mrp) {
	return mrp_to_dcm(dcm_to_quat(mrp));
}


arma::mat tilde(const arma::vec & x) {
	arma::mat vec_tilde = {
		{0, - x(2), x(1)},
		{x(2), 0, -x(0)},
		{ - x(1), x(0), 0}
	};

	return vec_tilde;
}

arma::mat M1(const double angle) {
	/**
	Returns the matrix of the elemental M1 rotation
	@param euler_angle Euler angle
	@return m1 Elemental rotation matrix
	*/
	arma::mat M = {
		{ 1, 0, 0},
		{0, cos(angle), sin(angle)},
		{0, - sin(angle), cos(angle)}
	};
	return M;
}


arma::mat M2(const double angle) {
	/**
	Returns the matrix of the elemental M2 rotation
	@param euler_angle Euler angle
	@return m2 Elemental rotation matrix
	*/
	arma::mat M = {
		{ cos(angle), 0, - sin(angle)},
		{0, 1, 0},
		{sin(angle), 0, cos(angle)}
	};
	return M;
}

arma::mat M3(const double angle) {
	/**
	Returns the matrix of the elemental M3 rotation
	@param euler_angle Euler angle
	@return m3 Elemental rotation matrix
	*/
	arma::mat M = {
		{cos(angle), sin(angle), 0},
		{ - sin(angle), cos(angle), 0},
		{0, 0, 1}
	};
	return M;
}

arma::vec dmrpdt(double t, arma::vec & attitude_set ) {
	arma::vec mrp = attitude_set.rows(0, 2);
	arma::vec omega = attitude_set.rows(3, 5);
	arma::mat I = arma::eye<arma::mat>(3, 3);

	return 0.25 * ( (1 - arma::dot(mrp, mrp)) * I + 2 * tilde(mrp) + 2 * mrp * mrp.t()) * omega;
}

arma::vec rot_dyn(double t, arma::vec & attitude_set, arma::mat & inertia) {
	arma::vec omega = attitude_set.rows(3, 5);
	arma::vec omega_dot = - arma::inv_sympd(inertia) * tilde(omega) * inertia * omega;
	return omega_dot;

}

arma::vec shadow_mrp(arma::vec mrp) {
	if (arma::norm(mrp) > 1) {
		return - mrp / arma::dot(mrp, mrp);
	}
	else {
		return mrp;
	}
}



// TO IMPLEMENT
// std::pair<double, arma::vec > quat_to_prv(const arma::vec & Q) {
// 	std::pair<double, arma::vec > principle_rotation_vector;


// 	return principle_rotation_vector;
// }

arma::vec quat_to_mrp(const arma::vec & Q , const bool short_rot) {
	arma::vec mrp = {Q(1), Q(2), Q(3)};
	mrp = mrp / ( 1 + Q(0));

	if (short_rot == true) {
		if (arma::norm(mrp) > 1) {
			mrp = - mrp / pow(arma::norm(mrp), 2);
		}
	}

	return mrp;
}

arma::vec dcm_to_quat(const arma::mat & dcm) {

	arma::vec Q = {0, 0, 0, 0};
	arma::vec q_s = {
		1. / 4. * ( 1 + arma::trace(dcm) ),
		1. / 4. * (1 + 2 * dcm(0, 0) - arma::trace(dcm)),
		1. / 4. * (1 + 2 * dcm(1, 1) - arma::trace(dcm)),
		1. / 4. * (1 + 2 * dcm(2, 2) - arma::trace(dcm))
	};


	int max_coef_index = q_s.index_max();
	Q(max_coef_index) = std::sqrt(q_s.max());

	switch (max_coef_index) {
	case 0:
		Q(1) = 1. / 4. * (dcm(1, 2) - dcm(2, 1)) / Q(0);
		Q(2) = 1. / 4. * (dcm(2, 0) - dcm(0, 2)) / Q(0);
		Q(3) = 1. / 4. * (dcm(0, 1) - dcm(1, 0)) / Q(0);
		break;

	case 1:
		Q(0) = 1. / 4. * (dcm(1, 2) - dcm(2, 1)) / Q(1);
		Q(2) = 1. / 4. * (dcm(0, 1) + dcm(1, 0)) / Q(1);
		Q(3) = 1. / 4. * (dcm(2, 0) + dcm(0, 2)) / Q(1);
		break;

	case 2:
		Q(0) = 1. / 4. * (dcm(2, 0) - dcm(0, 2)) / Q(2);
		Q(1) = 1. / 4. * (dcm(0, 1) + dcm(1, 0)) / Q(2);
		Q(3) = 1. / 4. * (dcm(1, 2) + dcm(2, 1)) / Q(2);
		break;


	case 3:
		Q(0) = 1. / 4. * (dcm(0, 1) - dcm(1, 0)) / Q(3);
		Q(1) = 1. / 4. * (dcm(2, 0) + dcm(0, 2)) / Q(3);
		Q(2) = 1. / 4. * (dcm(1, 2) + dcm(2, 1)) / Q(3);
		break;

	}

	return Q;

}

std::pair<double, arma::vec > dcm_to_prv(const arma::mat & dcm) {
	double angle = std::acos(0.5 * (arma::trace(dcm) - 1));
	arma::vec axis = {
		dcm(1, 2) - dcm(2, 1),
		dcm(2, 0) - dcm(0, 2),
		dcm(0, 1) - dcm(1, 0)
	};
	axis = 1. / (2 * std::sin(angle)) * axis;
	std::pair<double, arma::vec > prv;
	prv.first = angle;
	prv.second = axis;
	return prv;
}

arma::vec dcm_to_mrp(const arma::mat & dcm, const bool short_rot) {
	return quat_to_mrp(dcm_to_quat(dcm), short_rot);

}

arma::vec dcm_to_euler321(const arma::mat & dcm) {
	arma::vec angles = {0, 0, 0};
	angles(0) = std::atan2(dcm(0, 1) , dcm(0, 0));
	angles(1) = - std::asin(dcm(0, 2));
	angles(2) = std::atan2(dcm(1, 2), dcm(2, 2));
	return angles;
}

arma::vec dcm_to_euler313(const arma::mat & dcm) {
	arma::vec angles = {0, 0, 0};
	angles(0) = std::atan2(dcm(2, 0) , - dcm(2, 1));
	angles(1) = std::acos(dcm(2, 2));
	angles(2) = std::atan2(dcm(0, 2), dcm(1, 2));

	return angles;
}

arma::vec mrp_to_euler321(const arma::vec & sigma) {
	return dcm_to_euler321(mrp_to_dcm(sigma));
}

arma::vec mrp_to_euler313(const arma::vec & sigma) {
	return dcm_to_euler313(mrp_to_dcm(sigma));
}


arma::vec euler321_to_mrp(const arma::vec & euler_angles) {
	return dcm_to_mrp(euler321_to_dcm(euler_angles));
}

arma::vec euler313_to_mrp(const arma::vec & euler_angles) {
	return dcm_to_mrp(euler313_to_dcm(euler_angles));
}

arma::vec euler321d_to_mrp(arma::vec & euler_angles) {
	return dcm_to_mrp(euler321d_to_dcm(euler_angles));
}

arma::vec euler313d_to_mrp(arma::vec & euler_angles) {
	return dcm_to_mrp(euler313d_to_dcm(euler_angles));
}

arma::vec dcm_to_euler321d(const arma::mat & dcm) {
	return 180. / arma::datum::pi * dcm_to_euler321(dcm);
}

arma::vec dcm_to_euler313d(const arma::mat & dcm) {
	return 180. / arma::datum::pi * dcm_to_euler313(dcm);
}

arma::vec mrp_to_euler313d(const arma::vec & sigma) {
	return dcm_to_euler313d(mrp_to_dcm(sigma));
}

arma::vec mrp_to_euler321d(const arma::vec & sigma) {
	return dcm_to_euler321d(mrp_to_dcm(sigma));
}

arma::vec prv_to_mrp(const arma::vec & prv) {
	return dcm_to_mrp(prv_to_dcm(prv));
}

arma::mat prv_to_dcm(const arma::vec & prv) {

	if (arma::norm(prv) > 0) {
		double Phi = arma::norm(prv);
		double Sigma = 1 - std::cos(Phi);
		arma::vec e = prv / Phi;

		arma::mat dcm = Sigma * e * e.t() + std::cos(Phi) * arma::eye<arma::mat>(3,3) - std::sin(Phi) * tilde(e);


		return dcm;
	}

	else {
		return arma::eye<arma::mat>(3,3);
	}


}







