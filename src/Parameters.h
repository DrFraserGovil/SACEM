#pragma once
#include <random>
#include <assert.h>
template<typename T> class Parameter
{
	public:
		Parameter()
		{
			Value = 0;
		}
		Parameter(T value)
		{
			Value = value;
		}	
		T Value;
				
};

template<typename T> class BoundedParameter : public Parameter<T>
{
	public:
		BoundedParameter()
		{
			this->Value = 0;
			MinValue = 0;
			MaxValue = 0;
		}
	
		BoundedParameter(T value, T min, T max)
		{
			this->Value = value; //the "this" is needed to force the compiler to recognise Value as belonging to the parent class
			MinValue = min;
			MaxValue = max;
		}
		BoundedParameter(T min, T max)
		{
			this->Value = min;
			MinValue = min;
			MaxValue = max;
		}
		T MinValue;
		T MaxValue;
	
};

template<typename T> class IterableParameter : public BoundedParameter<T>
{
	public:
		IterableParameter(){};
		IterableParameter(T value, T min, T max, int steps)
		{
			this->Value = value;
			this->MinValue = min;
			this->MaxValue = max;
			NSteps = steps;
			
			CurrentIndex = 0;
			if (NSteps > 1)
			{
				bruch = (max - min)/(steps - 1);
			}
			else
			{
				bruch = 0;
				
				//if Nsteps < 0, just force everything into the midpoint of the provided values.
				this->MinValue = (max + min)/2;
				this->MaxValue = this->MinValue;
				this->Value = this->MinValue;
			}
		}
		
		void IterateValue()
		{
			++CurrentIndex;
			if (CurrentIndex >= NSteps)
			{
				CurrentIndex = 0;
			}
			this->Value = this->MinValue + bruch * CurrentIndex; 
		}
		void IterateValue(int stepIndex)
		{
			this->Value = this->MinValue + bruch*stepIndex;
			CurrentIndex = stepIndex;
		}
		
		T getStepValue(int stepIndex)
		{
			return this->MinValue + bruch*stepIndex;
		} 
		
		
		int NSteps;
	
	private:
		double bruch;
		int CurrentIndex;
};



template<typename T> class RandomisableParameter : public BoundedParameter<T>
{
	public:
		RandomisableParameter()
		{
			Initialised = false;
		}
		RandomisableParameter(T value, T min, T max, std::mt19937  * engine)
		{
			this->Value = value;
			this->MinValue = min;
			this->MaxValue = max;
			Randomiser = std::uniform_real_distribution(0.0,1.0);
			Engine = engine;
			Initialised = true;
		}
		RandomisableParameter(T min, T max, std::mt19937 * engine)
		{
			this->Value = min;
			this->MinValue = min;
			this->MaxValue = max;
			Randomiser = std::uniform_real_distribution(0.0,1.0);
			Engine = engine;
			Initialised = true;
		}
		
		
		void Scramble()
		{
			bool init = Initialised;
			assert(init);
			double bruch = Randomiser(*Engine);
			this->Value = this->MinValue + bruch * (this->MaxValue - this->MinValue);
		}
		
		
	private:
	
		std::uniform_real_distribution<double> Randomiser;
		std::mt19937 * Engine;
		bool Initialised;
};




