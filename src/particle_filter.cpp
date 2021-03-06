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
using std::uniform_int_distribution;
using std::uniform_real_distribution;


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
  num_particles = 100;  // Experimenting with 50 particles as of now.

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
  // std::cout << "Particles initialized! :) \n";


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
  // std::cout << "Particles prediction step done! :) \n";

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
    double min_distance = std::numeric_limits<double>::max();    // For each observation setting very high value then setting minimum value in the following for loop
    for(LandmarkObs& PRED: predicted){
      double distance = dist(OBS.x, OBS.y, PRED.x, PRED.y);
      if(distance < min_distance){	
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

  // for each particle specified
  for (int i = 0; i < num_particles; i++) {

    // acquire the particle x, y coordinates
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    // create a vector to hold the map landmark locations predicted to be within sensor range of the particle
    vector<LandmarkObs> predictions_vector;

    // for each map landmark specified
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

      // get x,y coordinates for each landmark
      float landmark_x = map_landmarks.landmark_list[j].x_f;
      float landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      
      // Consider landmarks within sensor range of the particle
      if (fabs(landmark_x - particle_x) <= sensor_range && fabs(landmark_y - particle_y) <= sensor_range) {

        // add prediction to vector
        predictions_vector.push_back(LandmarkObs{ landmark_id, landmark_x, landmark_y });
      }
    }

    // create and populate a copy of the list of observations transformed from vehicle coordinates to map coordinates
    //formulas are referred from Landmarks Quiz of Implementation of particle filter lesson.
    vector<LandmarkObs> transformed_obs;
    for (unsigned int j = 0; j < observations.size(); j++) {
      double transformed_x = particle_x + cos(particle_theta)*observations[j].x - sin(particle_theta)*observations[j].y;
      double transformed_y = particle_y + sin(particle_theta)*observations[j].x + cos(particle_theta)*observations[j].y;
      transformed_obs.push_back(LandmarkObs{ observations[j].id, transformed_x, transformed_y });
    }

    // perform dataAssociation for the predictions and transformed observations on i th particle
    dataAssociation(predictions_vector, transformed_obs);

    // reinitialize weight to 1.0
    particles[i].weight = 1.0;

    for (unsigned int j = 0; j < transformed_obs.size(); j++) {
      
      // observation and associated prediction coordinates
      double obs_x, obs_y, pred_x, pred_y;
      obs_x = transformed_obs[j].x;
      obs_y = transformed_obs[j].y;

      int associated_prediction = transformed_obs[j].id;

      // get the x,y coordinates of the prediction associated with the current observation
      for (unsigned int k = 0; k < predictions_vector.size(); k++) {
        if (predictions_vector[k].id == associated_prediction) {
          pred_x = predictions_vector[k].x;
          pred_y = predictions_vector[k].y;
        }
      }

      // calculate weight for this observation with multivariate Gaussian
      // formulas are referred from Particle Weights Solution section of Implementation of particle filter lesson.
      double std_lm_x = std_landmark[0]; // Landmark measurement uncertainty for x
      double std_lm_y = std_landmark[1]; // Landmark measurement uncertainty for y
      double obs_w = ( 1/(2*M_PI*std_lm_x*std_lm_y)) * exp( -( pow(pred_x-obs_x,2)/(2*pow(std_lm_x, 2)) + (pow(pred_y-obs_y,2)/(2*pow(std_lm_y, 2))) ) );

      // Resultant weight is product of this obersvation weight with total observations weight
      particles[i].weight *= obs_w;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<Particle> resampled_particles;

  // acquire all particles current weight and storing it in a vector
  vector<double> all_weights;
  for (int i = 0; i < num_particles; i++) {
    all_weights.push_back(particles[i].weight);
  }

  // Generate random starting index for resampling
  uniform_int_distribution<int> dist_uniform_int(0, num_particles-1);
  int index = dist_uniform_int(gen);

  // Get maximum weight from the all_weights vector
  double MW = *max_element(all_weights.begin(), all_weights.end());

  // uniform random distribution [0.0, Maximum_weight]
  uniform_real_distribution<double> dist_uniform_real(0.0, MW);

  // Initialize beta : Referred from Resampling Wheel section of particle filter lesson
  double beta = 0.0;

  // Referred from Resampling Wheel section of particle filter lesson
  for (int i = 0; i < num_particles; i++) {
    beta = beta + dist_uniform_real(gen) * 2.0; // multiplied by two for adding distribution in both directions
    while ( all_weights[index] < beta ) {
      beta = beta - all_weights[index];
      index = (index + 1) % num_particles; // % num_particles to avoid out of bound and cyclic world assumption :)
    }
    resampled_particles.push_back(particles[index]);
  }

  particles = resampled_particles; // Set resampled particles for the next step

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