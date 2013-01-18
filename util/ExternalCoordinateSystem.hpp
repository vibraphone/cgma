
#ifndef EXTERNAL_COORDINATE_SYSTEM_HPP
#define EXTERNAL_COORDINATE_SYSTEM_HPP

class ExternalCoordinateSystem
{
  public:

  virtual ~ExternalCoordinateSystem() {}

	//given s,n return x,y,z
	virtual bool GetXYZ( double i, double j, double k,
   		    		 double &x, double &y, double &z ) = 0;
   	
  //given x,y,z , return an s,n
  virtual bool GetIJK( double x, double y, double z,
	      			double &i, double &j, double &k ) = 0;	  
};

#endif

