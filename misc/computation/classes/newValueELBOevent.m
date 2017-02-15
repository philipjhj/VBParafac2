classdef (ConstructOnLoad) newValueELBOevent < event.EventData
   properties
      newValue
   end
   
   methods
      function data = newValueELBOevent(newValue)
         data.newValue = newValue;
      end
   end
end