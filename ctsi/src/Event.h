#ifndef EVENT_H
#define EVENT_H

#include <list>

using std::list;

class Event{
	public:
		Event();
		~Event();
		void          setEventID(int eventIDArg);
		void          addHit(float timeArg, float energyArg, float xArg, float yArg, float zArg, int detectorIDArg);
		void          flushEvent();
		int           getEventID();
		list<float>   getTime();
		void          setTime(list<float> timeIn);
		list<float>   getEnergy();
		void          setEnergy(list<float> energyIn);
		list<float>   getX();
		void          setX(list<float> xIn); 
		list<float>   getY();
		void          setY(list<float> yIn);
		list<float>   getZ();
		void          setZ(list<float> zIn);
		list<int>     getDetectorID();
		void          setDetectorID(list<int> detectorIDIn);
		list<int>     getAnode();
		void          setAnode(list<int> anodeIn);
		list<int>     getCathode();
		void          setCathode(list<int> cathodeIn);
		int           getNumIntrxns();

	protected:
		int           eventID;
		list<float>   time;
		list<float>   energy;
		list<float>   x;
		list<float>   y;
		list<float>   z;
		list<int>     detectorID;
		list<int>     anode;
		list<int>     cathode;
};

#endif //EVENT_H
