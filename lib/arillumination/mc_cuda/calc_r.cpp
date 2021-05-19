/* File: calc_r.cpp
 * Author: Urs Hofmann
 * Mail: hofmannu@biomed.ee.ethz.ch
 * Date: 02.07.2019
 * Description: Calculates the reflectance for different angles.
 */

#include <math.h>

void calc_r(
	float* RVector, // output reflectance vector
	const unsigned int n, // number of discretization steps for the angle
	const float dalpha, // step size of angle [rad]
	const float nWater, // refractive index of water
	const float nTissue // refractive index of tissue
){
	
	float alpha; // angle of incident beam
	float beta; // angle of outgoing beam
	float rs; // reflection coefficient for parallel pol
	float rp; // reflection coefficient for vertically pol
	float* RVector_cp = RVector;
	float minR = 1000;
	float maxR = 0;
	float crit_angle = asin(nWater / nTissue);
	
	for (unsigned int ir = 0; ir < n; ir++){
	
		// calculate incoming angle for current step
		alpha = dalpha * ir;
		
		// use snells law to calculate outgoing angle
		beta = asin(nTissue / nWater * sin(alpha));

		rp = tan(alpha - beta) / tan(alpha + beta);
		rs = sin(alpha - beta) / sin(alpha + beta); 	
		if (alpha == 0)
			*RVector_cp = pow(((nTissue - nWater) / (nTissue + nWater)), 2);
		else if (alpha > crit_angle)
			*RVector_cp = 1; // total reflection
		else
			*RVector_cp = 0.5 * (pow(rp, 2) + pow(rs, 2));

		// printf("%2.2f %2.2f %2.2f\n", alpha, beta, *RVector_cp);
		// check for min and max value
		if ((*RVector_cp) > maxR)
			maxR = *RVector_cp;

		if ((*RVector_cp) < minR)
			minR = *RVector_cp;

		RVector_cp++; // shift pointer to next position
	}

	// printf("R ranges from %2.6f to %2.6f (dalpha = %f) \n", minR, maxR, dalpha);	

}
