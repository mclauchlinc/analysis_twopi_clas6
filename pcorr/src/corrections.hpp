#ifndef CORRECTIONS_HPP
#define CORRECTIONS_HPP

#include "constants.hpp"
#include "functions.hpp"
#include <math.h>
#include <stdio.h>

namespace corr{
	//{e16/e1f} {Sector} {A,B,C,D,E} {alpha,beta,gamma}
	static const float _angle_e_par[2][6][5][3] = 		{	{	{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}}},
											{	{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
												{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}}}};
	//{e16,e1f} {sector} {A,B,C,D} {alpha, beta, gamma}
	static const float _p_e_par[2][6][4][3] = 	{	{	{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}}},
										{	{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}},
											{	{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3},{0.1,0.2,0.3}}}};

	float p_corr_e(float p_e_, float theta_e_, float phi_e_, int run_, bool centered_ = false, int sector_=0);
	float theta_e_corr(float theta_e_, float phi_e_, int run_, bool centered_ = false, int sector_=0);
}


#endif
/*
{}
{x,x}//{e16,e1f}
{{x,x,x,x,x,x},{x,x,x,x,x,x}}//{e16,e1f} {sector}
{{{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x}},{{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x},{x,x,x,x,x}}}//{e16,e1f} {sector} {A,B,C,D,E}
{{{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}}},{{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},{{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}}}}//{e16,e1f} {sector} {A,B,C,D,E} {alpha,beta,gamma}

{	{	{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}}},
	{	{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}},
		{	{x,x,x},{x,x,x},{x,x,x},{x,x,x},{x,x,x}}}}//{e16,e1f} {sector} {A,B,C,D,E} {alpha,beta,gamma}
*/
