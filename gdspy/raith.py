import struct
import sys
import numpy
from gdspy.gdsiiformat import (
    _eight_byte_real,
)
#ported from https://github.com/ahryciw/Raith_GDSII/blob/main/src/Raith_library.m

class RaithEllipse():
    def __init__(self,center,radMajMin,linewidth=None,number_of_points=20,angle=None,layer=0,datatype=0):
        if number_of_points>1024:
            raise ValueError("Does not support more than 1024 vertices")
        self.number_of_points = number_of_points
        self.layer            = layer
        self.datatype         = datatype
        self.center           = center
        self.radMajMin        = radMajMin #radius major minor axis, use array with same value twice for circle
        if self.radMajMin[0] == self.radMajMin[1]:
            self.angle            = None
        else:
            self.angle            = angle #only makes sense for ellipse
        self.linewidth        = linewidth #in um, when =0 then fill circle
    
    def to_gds(self, outfile, multiplier):
        bytesString = bytes()
        flags = 0
        if self.radMajMin[0] != self.radMajMin[1]:
            flags += 1
        bytesString += struct.pack(">HH HHh HHh",
                4,0x5600,
                6,0x0D02, self.layer,          #layer
                6,0x0E02, self.datatype,       #dose factor * 1000
        )
        if self.linewidth is not None:
            bytesString += struct.pack(">HHi",
                8,0x0F03, int(self.linewidth*1000), #linewidth (convert to nm)
            )
        else:
            flags += 2 #filled
        if self.angle is not None:
            bytesString += struct.pack(">HH",12,0x1c05)
            bytesString += _eight_byte_real(self.angle*1.0)#in deg relative to u axis          
        bytesString += struct.pack(">HH8i HH",
                36,0x1003,int(self.center[0]*1000), int(self.center[1]*1000), 
                          int(self.radMajMin[0]*1000), int(self.radMajMin[1]*1000),
                          0, 0, 
                          self.number_of_points, flags,
                4,0x1100 #end element
        )
        outfile.write(bytesString)
            
            
        
class RaithArc():
    def __init__(self,center,radMajMin,startEndAngle, linewidth=None,number_of_points=20,angle=None,layer=0,datatype=0):
        if number_of_points>1024:
            raise ValueError("Does not support more than 1024 vertices")
        self.number_of_points = number_of_points
        self.layer            = layer
        self.datatype         = datatype
        self.center           = center
        self.radMajMin        = radMajMin #radius major minor axis, use array with same value twice for circle
        self.startEndAngle    = startEndAngle #angle from u axis, counterclockwise degree
        if self.radMajMin[0] == self.radMajMin[1]:
            self.angle            = None
        else:
            self.angle            = angle #only makes sense for ellipse
        self.linewidth        = linewidth #in um, when =0 then fill circle
    
    def to_gds(self, outfile, multiplier):
        bytesString = bytes()
        flags = 4
        if self.radMajMin[0] != self.radMajMin[1]:
            flags += 1
        bytesString += struct.pack(">HH HHh HHh",
                4,0x5600,
                6,0x0D02, self.layer,          #layer
                6,0x0E02, self.datatype,       #dose factor * 1000
        )
        if self.linewidth is not None:
            bytesString += struct.pack(">HHi",
                8,0x0F03, int(self.linewidth*1000), #linewidth (convert to nm)
            )
        else:
            flags += 2 #filled
        if self.angle is not None:
            bytesString += struct.pack(">HH",12,0x1c05)
            bytesString += _eight_byte_real(self.angle*1.0)#in deg relative to u axis          
        bytesString += struct.pack(">HH8i HH",
                36,0x1003,int(self.center[0]*1000), int(self.center[1]*1000), 
                          int(self.radMajMin[0]*1000), int(self.radMajMin[1]*1000),
                          int(785398/45*self.startEndAngle[0]), int(785398/45*self.startEndAngle[1]), 
                          self.number_of_points, flags,
                4,0x1100 #end element
        )
        outfile.write(bytesString)


class FbmsCircle():
    #zero linewidth means single pixel line
    def __init__(self, centerPoint, radiusUm, linewidthUm = 0, layer=0, datatype=0):
        self._polygon_dict = None
        self.widths = numpy.array([[radiusUm, linewidthUm]])
        self.points    = numpy.array([centerPoint])
        self.layers    = [layer]
        self.datatypes = [datatype]
        self.gdsii_path = True
        self.properties = {}
    def to_gds(self, outfile, multiplier):
        bytesStart = struct.pack(">4Hh2Hh2Hl",
                                4,0x5800,
                                6,0x0D02,self.layers[0],
                                6,0x0E02,self.datatypes[0],
                                8,0x0F03,int(round(self.widths[0, 1] * multiplier)))
        outfile.write(bytesStart)
        points = numpy.round(self.points * multiplier)
        uvData = numpy.zeros(8, dtype=">i4")
        uvData[5] = points[:,0]
        uvData[6] = points[:,1]
        uvData[7] = numpy.round(self.widths[0, 0] * multiplier)
        outfile.write(struct.pack(">2H", 4 + 4 * uvData.size, 0x1003))
        outfile.write(uvData.tobytes())
        outfile.write(struct.pack(">2H", 4, 0x1100))

            
class FbmsPath():
    #segmentsCurvature in um, positive value = right curve
    def __init__(self, points, width, segmentsCurveDistCenter=0, tolerance=0.01, precision=1e-3, max_points=199, layer=0, datatype=0):
        self._polygon_dict = None
        self.n = 1
        self.points = numpy.array(points)
       
        self.segmentsCurveDistCenter = numpy.array(segmentsCurveDistCenter)
        if self.segmentsCurveDistCenter.size == 1:
            self.segmentsCurveDistCenter = numpy.repeat(self.segmentsCurveDistCenter, self.points.shape[0]-1, axis=0)
        if self.segmentsCurveDistCenter.size != (self.points.shape[0] - 1):
            raise ValueError("segmentsCurveDistCenter is underdefined")
        self.widths = numpy.array([[width] + list(self.segmentsCurveDistCenter)])
        
        self.layers = [layer]
        self.datatypes = [datatype]
        self.tolerance = tolerance
        self.precision = precision
        self.max_points = max_points
        self.gdsii_path = True
        self.properties = {}
    
    #information from https://github.com/ahryciw/Raith_GDSII/blob/cbe3929786d646d09ea9449ec1f666a4afdadf91/src/Raith_library.m#L793-L807
    def to_gds(self, outfile, multiplier):
        if len(self.points) < 2:
            return
        if self.n != 1:
            return
        bytesStart = struct.pack(">4Hh2Hh2Hl",
                                4,0x5800,
                                6,0x0D02,self.layers[0],
                                6,0x0E02,self.datatypes[0],
                                8,0x0F03,int(round(self.widths[0, 0] * multiplier)))
        outfile.write(bytesStart)
        
        segmentType = (self.widths[0, 1:] != 0).astype("int") + 1
        
        points = numpy.round(self.points * multiplier)
        uvData = numpy.zeros(points.size*2+4, dtype=">i4")
        
        uvData[8::4] = segmentType
        uvData[5::4] = points[:,0]
        uvData[6::4] = points[:,1]
        uvData[11::4] = numpy.round(self.widths[0, 1:] * multiplier)
        outfile.write(struct.pack(">2H", 4 + 4 * uvData.size, 0x1003))
        outfile.write(uvData.tobytes())
        outfile.write(struct.pack(">2H", 4, 0x1100))
