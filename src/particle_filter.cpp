/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	is_initialized = false;
	std::default_random_engine gen;

	num_particles = 100;
	
	// Set standard deviations for x, y, and psi.
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	
	// Create normal distributions for x, y and psi.
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	weights.resize(num_particles);
	particles.resize(num_particles);

	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;
		
		// sample_x = dist_x(gen);
		// where "gen" is the random engine initialized earlier (line 18).
        particles[i].id = i;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
		particles[i].weight = 1/num_particles;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;

	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	std::normal_distribution<double> dist_x(0, std_x);
	std::normal_distribution<double> dist_y(0, std_y);
	std::normal_distribution<double> dist_theta(0, std_theta);

	double theta_next;
	// Create normal distributions for x, y and psi.
	
	for (int i = 0; i < num_particles; ++i) {
		if (std::fabs(yaw_rate) < 0.001){
			particles[i].x += velocity * delta_t * std::cos(particles[i].theta);
            particles[i].y += velocity * delta_t * std::sin(particles[i].theta);
		}
		else{	
			theta_next = particles[i].theta + yaw_rate * delta_t;
			particles[i].x += (velocity/yaw_rate) * ( std::sin(theta_next) - std::sin(particles[i].theta));
        	particles[i].y += (velocity/yaw_rate) * (-std::cos(theta_next) + std::cos(particles[i].theta));
        	particles[i].theta = theta_next;
        }
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	
	//Using method described by driveWell (https://discussions.udacity.com/t/implementing-data-association/243745/2):
	//	for each particle:
    //	  for each observation:
    //    	  transform_observation_to_map
    //    	  for each landmark:
    //        	  calculate euclidean distance and associate to particle the landmark id with the min distance
	double std_x = std_landmark[0];
    double std_y = std_landmark[1];

	double total_multivar_prob;
	double tmp_prob;

	double x_observation;
	double y_observation;

	double tmp_landmark_distance;
	int min_landmark_id;
	double min_landmark_distance;

	bool landmark_in_sensor_range;

	double d_x_squre;
	double d_y_squre;

	double total_weight = 0;

	for (int i = 0; i < num_particles; ++i) {
		total_multivar_prob = 1.0;

		for (int j = 0; j < observations.size(); ++j) {
			//transform_observation_to_map
			x_observation = particles[i].x + observations[j].x * std::cos(particles[i].theta) - observations[j].y * std::sin(particles[i].theta);
            y_observation = particles[i].y + observations[j].x * std::sin(particles[i].theta) + observations[j].y * std::cos(particles[i].theta);

            min_landmark_distance = sensor_range;
            landmark_in_sensor_range = false;

            //finding the associated landmark
			for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) 
			{
				//calculate euclidean distance and associate to particle the landmark id with the min distance

				tmp_landmark_distance = dist(x_observation, y_observation, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
				
				if (tmp_landmark_distance <= min_landmark_distance){
					landmark_in_sensor_range = true;
					min_landmark_distance = tmp_landmark_distance;
					min_landmark_id = k;
				}
			}
			if (not landmark_in_sensor_range){
				continue;
			}
			//associated landmark found
			//calculating multivariate_gaussian_prob
			d_x_squre = pow((x_observation - map_landmarks.landmark_list[min_landmark_id].x_f), 2);
			d_y_squre = pow((y_observation - map_landmarks.landmark_list[min_landmark_id].y_f), 2);
			tmp_prob = (1/(2*M_PI*std_x*std_y)) * exp(-( d_x_squre/(2*pow(std_x,2)) + d_y_squre/(2*pow(std_y,2)) )) ;
			
			total_multivar_prob *= tmp_prob;
		}

		//updating un-normalized weight
		if (total_multivar_prob == 1) {
			particles[i].weight = 0;
		}
		else {
			particles[i].weight = total_multivar_prob;
		}
		total_weight += particles[i].weight;
		
	}

	//normalizing wieght and updating weights[] 
	for (int i = 0; i < num_particles; ++i) {
		particles[i].weight = particles[i].weight / total_weight;
        weights[i] = particles[i].weight;
    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    std::vector<Particle> particles_curr(particles);
    std::default_random_engine gen;
    std::discrete_distribution<std::size_t> dis(weights.begin(), weights.end());

    for (int i = 0; i < num_particles; ++i) {
        particles[i] = particles_curr[dis(gen)];
    }
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
