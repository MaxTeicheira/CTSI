#include "Event.h"

using namespace std;

Event::Event() : eventID(0)
{
}

Event::~Event()
{
}

void Event::setEventID(int eventIDArg)
{
	eventID = eventIDArg;
}

void Event::addHit(float timeArg, float energyArg, float xArg, float yArg, float zArg, int detectorIDArg)
{
	time.push_back(timeArg);
	energy.push_back(energyArg);
	x.push_back(xArg);
	y.push_back(yArg);
	z.push_back(zArg);
	detectorID.push_back(detectorIDArg);
}

void Event::flushEvent()
{
	eventID = 0;
	time.clear();
	energy.clear();
	x.clear();
	y.clear();
	z.clear();
	detectorID.clear();
}

int Event::getEventID()
{
	return eventID;
}

list<float> Event::getTime()
{
	return time;
}

void Event::setTime(list<float> timeIn)
{
	time = timeIn;
}

list<float> Event::getEnergy()
{
	return energy;
}

void Event::setEnergy(list<float> energyIn)
{
	energy = energyIn;
}

list<float> Event::getX()
{
	return x;
}

void Event::setX(list<float> xIn)
{
	x = xIn;
}

list<float> Event::getY()
{
	return y;
}

void Event::setY(list<float> yIn)
{
	y = yIn;
}

list<float> Event::getZ()
{
	return z;
}

void Event::setZ(list<float> zIn)
{
	z = zIn;
}

list<int> Event::getDetectorID()
{
	return detectorID;
}

void Event::setDetectorID(list<int> detectorIDIn)
{
	detectorID = detectorIDIn;
}

list<int> Event::getAnode()
{
	return anode;
}

void Event::setAnode(list<int> anodeIn)
{
	anode = anodeIn;
}

list<int> Event::getCathode()
{
	return cathode;
}

void Event::setCathode(list<int> cathodeIn)
{
	cathode = cathodeIn;
}

int Event::getNumIntrxns()
{
	return (int) energy.size();
}
