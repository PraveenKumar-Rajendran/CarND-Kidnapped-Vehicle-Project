/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>  // // Need this for sampling from distributions
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
std::default_random_engine gen;


  /**
   * init Initializes particle filter by initializing particles to Gaussian
   *   distribution around first position and all the weights to 1.
   * @param x Initial x position [m] (simulated estimate from GPS)
   * @param y Initial y position [m]
   * @param theta Initial orientation [rad]
   * @param std[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // TODO: Set the number of particles
  num_particles = 50;  // Experimenting with 50 particles as of now.

  double std_x = std[0]; // standard deviation of x [m]
  double std_y = std[1]; // standard deviation of y [m]
  double std_theta = std[2]; // standard deviation of yaw [rad]

  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std_x);
  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_y(y, std_y);
  // This line creates a normal (Gaussian) distribution for theta
  normal_distribution<double> dist_theta(theta, std_theta);
  
  //Took Reference from Program Gaussian Sampling from Implementation of Particle filter lesson
  //sampling from the above normal distributions for each particle.
  for (int i = 0; i < num_particles; ++i) {
    Particle single_particle;

    single_particle.id = i;
    single_particle.x = dist_x(gen);
    single_particle.y = dist_y(gen);
    single_particle.theta = dist_theta(gen);
    single_particle.weight = 1.0;  //Randomly initializing it to 1

    particles.push_back(single_particle); //Particle Element detail is added to the vector at the end.
  }

  is_initialized = true; //Making sure that init function is never run again for the next steps.
  std::cout << "Particles initialized! :) \n";


}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // Random Gaussian noise to contribute to sensor noise.
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  //Formula is referred from Yaw Rate and Velocity section of Motion model lesson of the nanodegree program
  // Add measurements to each particle
  for (int i = 0; i < num_particles; i++) {
    if (fabs(yaw_rate) >= 0.00001) {
      particles[i].x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta = particles[i].theta + yaw_rate * delta_t;
    }
    else {
      particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
    }
  // Add random Gaussian noise to each particle
  particles[i].x = particles[i].x + dist_x(gen);
  particles[i].y = particles[i].y + dist_y(gen);
  particles[i].theta = particles[i].theta + dist_theta(gen);
  std::cout << "Particles prediction step done! :) \n";

}

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  // Data Association: using Nearest Neighbor technique [ NOTE: works well in most cases however it is not efficient when particle density and Map density is high. ]
  for(LandmarkObs& OBS: observations){
    double min_distance = numeric_limits<double>::max();    // For each observation setting very high value then setting minimum value in the following for loop
    for(LandmarkObs& PRED: predicted){
      double distance = dist(OBS.x, OBS.y, PRED.x, PRED.y);
      if(d < min_distance){	
	min_distance = distance;
	OBS.id = PRED.id;  //Associating each observation to the closest landmark(nearest neighbour)
  // NOTE : //complexity is o(no of observations * no of predictions );
      }
    }    
  }  


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::discrete_distribution<> distribution_weights(weights.begin(), weights.end());
  vector<Particle> resampled_Particles;
  
  // NOTE: particles are sampled in such a way that it is proportional to its weight. ie, More the weight more the probability that it will get sampled repeatedly.
  //Loop over number of particles to create new ones with their probability of their weight
  for(int i = 0; i < num_particles; i++) {
    int selected_index_weight_prob = distribution_weights(gen);
    resampled_Particles.push_back(particles[selected_index_weight_prob]);
  }

  // Clear vectors
  particles.clear();
  weights.clear();

  particles = resampled_Particles; //Assign new set of particles with the resampled particles.
  
  std::cout<<"Finished Resampling! :) \n";

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}