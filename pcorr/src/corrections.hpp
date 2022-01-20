#ifndef CORRECTIONS_HPP
#define CORRECTIONS_HPP

#include "constants.hpp"
#include "functions.hpp"

namespace corr{
	//{e16/e1f} {Sector} {A,B,C,D,E} {alpha,beta,gamma}
	float _angle_e_par = 	{	{	{	{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
							{	{	{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}}}}};
	//{e16,e1f} {sector} {A,B,C,D} {alpha, beta, gamma}
	float _p_e_par = 	{	{	{	{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
							{	{	{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}},
									{	{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma},{alpha,beta,gamma}}}}};

	float p_corr_e(float p_e_, float theta_e_, float phi_e_, bool centered_ = false, int sector_=0);
	float theta_e_corr(float theta_e_, float phi_e_, bool centered_ = false, int sector_=0);
}


#endif