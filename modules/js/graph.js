/**
 * graph.js
 * graph.js is a library of functions necessary
 * to draw a 3d graph in Javascript/HTML
 */

/**
 * draws a 3d graph
 */
var displayGraphFunc = function() 
	{ 

		var f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
		graphData.set("expr", 
		    function (emit, x, y, i, j, time) 
			{ 
          var z = f(x,y)
      
			    emit( x, f(x,y), y );
			}
		);
		
		//autofit the scale of the graph so it fits to the screen nicely
		if (zAutofit)
		{
			var xStep = 1;// (xMax - xMin) / 256;
			var yStep = 1;//(yMax - yMin) / 256;
			var zSmallest = f(xMin, yMin);
			var zBiggest  = f(xMax, yMax);
      
      //infinity/NaN check
      if(!isFinite(zSmallest) || isNaN(zSmallest))
      {
        zSmallest = (xMin + yMin)/2;
      }
      
      if(!isFinite(zBiggest) || isNaN(zBiggest))
      {
        zSmallest = (xMax + yMax)/2;
      }

			//for every possible x & combination, find the smallest & biggest z values
			//infinity check is a must, otherwise when the equation results in infinity (eg divide by zero) then the graph cannot be plotted
			//when the z value is infinite, we can assign our own value to scale the graph nicely (this is optional)

			for (var x = xMin; x <= xMax; x += xStep)
			{
				for (var y = yMin; y <= yMax; y += yStep)
				{
					var z = f(x,y);
          
          //infinity/NaN check
          if(!isFinite(z) || isNaN(z))
          {
            continue;
          }
          
					if (z <= zSmallest){
						zSmallest = z;
					}

					if (z >= zBiggest) {
						zBiggest = z;
					}
				}
			}

			// change the global zMin var value to follow the smallest possible z val and
			// global zMax var value to follow the biggest possible z val
			zMin = zSmallest;
			zMax = zBiggest;
		}
		
	
		view.set("range", [[xMin, xMax], [zMin,zMax], [yMin, yMax]]); 
		
		//set the graph color to greenish hue
		graphColors.set("expr", 
				function (emit, x, y, i, j, time) 
				{ 
					var z = f(x,y);
          
					//var percent = (z - 1.2 * zMin) / (zMax - 1.2 * zMin);
          var percent = (z - zMin)/(zMax - zMin);
					var color = new THREE.Color( 0xffffff );
					color.setHSL( 1-percent, 1, 0.5);
          
          //hide emissions approaching infinity
          var alpha = 0.8;
          
          /*if(percent > 1.5 || percent < -0.5)
          {
            alpha = 0;
          }*/
          
					emit( color.r, color.g, color.b, alpha ); //last one is the alpha (transparency) setting
          //emit( 0.5, 0.8, percent, 0.4 );
          
				}
			);
	}

/**
 * updates the graph's tangents at specified points
 */	
var displayTangentsFunc = function ()

{

  var f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
  
  var xtFunc = Parser.parse( xtFunctionText ).toJSFunction( ['t'] );
  var ytFunc = Parser.parse( ytFunctionText ).toJSFunction( ['t'] );
  
  //set the point to reflect x and y functions of t
  xPoint = xtFunc(t); //x0 = x(t0)
  yPoint = ytFunc(t); //y0 = y(t0)
  zPoint = f(xPoint, yPoint); //f(x0,y0) = f(x(t0),y(t0))
  originPointData.set("data", [ [xPoint, zPoint, yPoint] ] );
  
  // set the point label data
  var pointXDisp = Math.round(xPoint * 100) / 100; //rounded to 2 decimals for better view
  var pointYDisp = Math.round(yPoint * 100) / 100;
  pointDataText.set( "data", ["(x(t0),y(t0)) = (" + pointXDisp + "," + pointYDisp + ")"] );
  
  
  // calculate the derivatives:
  
  //calculate the derivative of x with respect to t
  var dx_dt = math.derivative(xtFunctionText, 't').evaluate({t: t}) ; 
  
  //calculate the derivative of y with respect to t
  var dy_dt = math.derivative(ytFunctionText, 't').evaluate({t: t}) ;
  
  //calculate the derivative of f(x,y) with respect to x
  var dz_dx = math.derivative(functionText, 'x').evaluate({x: xPoint, y: yPoint}) ; 
  
  //calculate the derivative of f(x,y) with respect to y
  var dz_dy = math.derivative(functionText, 'y').evaluate({x: xPoint, y: yPoint}) ; 
  
  // calculate df of f(x,y) aka dz/dt
  var df = dz_dx * dx_dt + dz_dy * dy_dt;
  
  
  // to plot the derivatives lines, we need to create vectors
  // all vectors should originate from the trace point
  var vectorTail = [xPoint, zPoint, yPoint];
  
  // amount to scale vector by
  var scale = 1
  
  // plotting the derivative of f(x,y) with respect to x
  // since y is held constant, slope vector is [1,0,dz_dx]
  // don't forget that mathbox flips the y and z axis so it becomes [1,dz_dx,0]
  var dxVectorHead = [xPoint + scale, zPoint + scale * dz_dx, yPoint];
  dxDerivData.set( "data", [dxVectorHead,vectorTail] );
  
  // plotting the derivative of f(x,y) with respect to y
  // since x is held constant, slope vector is [0,1,dz_dx]
  // don't forget that mathbox flips the y and z axis so it becomes [0,dz_dx,1]
  var dyVectorHead = [xPoint, zPoint + scale * dz_dy, yPoint + scale];
  dyDerivData.set( "data", [dyVectorHead,vectorTail] );
  
  // calculate tangent plane from x and y slopes (dz_dx and dz_dy)
  tangentPlaneData.set("expr", 
      function (emit, x, y, i, j, time) 
    { 
      var z = dz_dx*(x - xPoint) + dz_dy*(y - yPoint) + zPoint;
        emit( x, z, y );
    }
  );
  
  //plot df as vector
  var dfVectorHead = [xPoint + scale, zPoint + scale * df, yPoint + scale];
  dfVectorData.set( "data", [dfVectorHead, vectorTail] );
  
  
}

/**
 * Toggle the colors of the graph
 */
var toggleGraphColorFunc = function()
{
  var f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
  
  if(colorGraph)
  {
    //set the graph color to rainbow hue
		graphColors.set("expr", 
				function (emit, x, y, i, j, time) 
				{ 
					var z = f(x,y);
					
					var percent = (z - 1.2 * zMin) / (zMax - 1.2 * zMin);
					var color = new THREE.Color( 0xffffff );
					color.setHSL( 1-percent, 1, 0.5 );
          
          //hide emissions approaching infinity
          var alpha = 0.4;
          
          /*if(percent > 1.5 || percent < -0.5)
          {
            alpha = 0;
          }*/
          
					emit( color.r, color.g, color.b, alpha ); //last one is the alpha (transparency) setting
				}
			);
  }
  else
  {
    //set the graph color to greenish hue
		graphColors.set("expr", 
				function (emit, x, y, i, j, time) 
				{ 
					var z = f(x,y);
					var percent = (z - zMin) / (zMax - zMin);
          
          //hide emissions approaching infinity
          var alpha = 0.8;
          
          /*if(percent > 1.5 || percent < -0.5)
          {
            alpha = 0;
          }*/
          
					emit( percent*0.6, percent, percent, alpha );
				}
			);	
  }
}

/**
 * updates the graph's partial derivatives at specified points
 */	
var displayDerivativesFunc = function ()
	
	{
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
		zPoint = f(xPoint, yPoint);
		tracePointData.set("data", [ [xPoint, zPoint, yPoint] ] );
		
		var xPointDisp = Math.round(xPoint * 100) / 100;
		var yPointDisp = Math.round(yPoint * 100) / 100;
		
		traceText.set("data", ["(" + xPointDisp + "," + yPointDisp + ")"]);

		xCurveData.set("expr", 
		    function (emit, x, i, t) 
			{ 
			    emit( x, f(x,yPoint), yPoint );
			}
		);

		yCurveData.set("expr", 
		    function (emit, y, j, t) 
			{ 
			    emit( xPoint, f(xPoint,y), y );
			}
		);
		
		
		// calculate partial derivatives
		
		// valid for all vectors
		var vectorTail = [xPoint, zPoint, yPoint];
		
		// amount to scale slope vector by
		var scale = 2
		
		//calculate the derivative of f(x,y) with respect to x
		var dz_dx = math.derivative(functionText, 'x').evaluate({x: xPoint, y: yPoint}) ; 
		
		//calculate the derivative of f(x,y) with respect to y
		var dz_dy = math.derivative(functionText, 'y').evaluate({x: xPoint, y: yPoint}) ; 
		
		// plotting the derivative of f(x,y) with respect to x
		// since y is held constant, slope vector is [1,0,dz_dx]
		// don't forget that mathbox flips the y and z axis so it becomes [1,dz_dx,0]
		var dxVectorHead = [xPoint + scale, zPoint + scale * dz_dx, yPoint];
		var dxVectorHead2 = [xPoint - scale, zPoint - scale * dz_dx, yPoint];
		dxData.set( "data", [dxVectorHead,vectorTail] );
		dxData2.set( "data", [dxVectorHead2,vectorTail] );
		
		// plotting the derivative of f(x,y) with respect to y
		// since x is held constant, slope vector is [0,1,dz_dx]
		// don't forget that mathbox flips the y and z axis so it becomes [0,dz_dx,1]
		var dyVectorHead = [xPoint, zPoint + scale * dz_dy, yPoint + scale];
		var dyVectorHead2 = [xPoint, zPoint - scale * dz_dy, yPoint - scale];
		dyData.set( "data", [dyVectorHead,vectorTail] );
		dyData2.set( "data", [dyVectorHead2,vectorTail] );
		
	}
	
 /**
  * updates the graph's directional derivatives at specified points
  */	
  var displayDirDerivativesFunc = function()
	{
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
		zPoint = f(xPoint, yPoint);
		tracePointData.set("data", [ [xPoint, zPoint, yPoint] ] );
    
		// calculate partial derivatives
		
		// valid for all vectors
		var vectorTail = [xPoint, zPoint, yPoint];
		
		// amount to scale slope vector by
		var scale = 1.5
		
		//calculate the derivative of f(x,y) with respect to x
		var dz_dx = math.derivative(functionText, 'x').evaluate({x: xPoint, y: yPoint}) ; 
		
		//calculate the derivative of f(x,y) with respect to y
		var dz_dy = math.derivative(functionText, 'y').evaluate({x: xPoint, y: yPoint}) ; 
		
		// plotting the derivative of f(x,y) with respect to x
		// since y is held constant, slope vector is [1,0,dz_dx]
		// don't forget that mathbox flips the y and z axis so it becomes [1,dz_dx,0]
		var dxVectorHead = [xPoint + scale, zPoint + scale * dz_dx, yPoint];

		dxData.set( "data", [dxVectorHead,vectorTail] );

		
		// plotting the derivative of f(x,y) with respect to y
		// since x is held constant, slope vector is [0,1,dz_dx]
		// don't forget that mathbox flips the y and z axis so it becomes [0,dz_dx,1]
		var dyVectorHead = [xPoint, zPoint + scale * dz_dy, yPoint + scale];

		dyData.set( "data", [dyVectorHead,vectorTail] );

		
		
		// calculate tangent plane function
		tangentPlaneData.set("expr", 
		    function (emit, x, y, i, j, t) 
			{ 
				var z = dz_dx*(x - xPoint) + dz_dy*(y - yPoint) + zPoint;
			    emit( x, z, y );
			}
		);
		
		
		// point (x0,y0) + unit vector (not projected on surface)
		// we ask the users to put in vector data, hence we need to calculate the unit vector for them
		var magnitude = Math.sqrt(Math.pow(vectorU,2)+Math.pow(vectorV,2));
		var unitVecU = vectorU / magnitude;
		var unitVecV = vectorV / magnitude;
		//unit vector has magnitude of 1, so zPoint is + 1
		var dirVectorHead = [xPoint + unitVecU, zPoint + 1, yPoint + unitVecV];
		dirVectorData.set( "data", [dirVectorHead,vectorTail] );
		
		
		// directional derivative (unit vector projected on surface/tangent plane)
		var dz = dz_dx * unitVecU + dz_dy * unitVecV;
		var gradientVectorHead = [xPoint + unitVecU * scale, zPoint + dz * scale, yPoint + unitVecV * scale];
		var gradientVectorHead2 = [xPoint - unitVecU * scale, zPoint - dz * scale, yPoint - unitVecV * scale];
		gradientData.set( "data", [gradientVectorHead, vectorTail] );
		gradientData2.set( "data", [gradientVectorHead2, vectorTail] );
		
    var dzDisp = Math.round(dz * 100) / 100;
    
    traceText.set("data", ["(" + dzDisp + ")"]);
		
	}
	
 /**
  * updates the graph's contours at specified points
  */
  var displayContourFunc = function ()
	{
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
		vertexArray = [];
		edgeArray = [];
		
		vertexArray = marchingSquares( xMin, xMax, yMin, yMax, f, c, 256 );
		
		for (var n = 0; n < vertexArray.length; n += 2)
		{
			edgeArray.push( vertexArray[n], vertexArray[n+1] );
		}
		edgeData.set( "width", edgeArray.length/2 );
		edgeData.set( "data", edgeArray );
		
		
		vertexArray2 = [];
		edgeArray2 = [];
		
		vertexArray2 = marchingSquares( xMin, xMax, yMin, yMax, f, c2, 256 );
		
		for (var n = 0; n < vertexArray2.length; n += 2)
		{
			edgeArray2.push( vertexArray2[n], vertexArray2[n+1] );
		}
		edgeData2.set( "width", edgeArray2.length/2 );
		edgeData2.set( "data", edgeArray2 );
		
		
		vertexArray3 = [];
		edgeArray3 = [];
		
		vertexArray3 = marchingSquares( cxMin, cxMax, cyMin, cyMax, f, c3, 256 );
		
		for (var n = 0; n < vertexArray3.length; n += 2)
		{
			edgeArray3.push( vertexArray3[n], vertexArray3[n+1] );
		}
		edgeData3.set( "width", edgeArray3.length/2 );
		edgeData3.set( "data", edgeArray3 );
		
	}
  
  /**
  * updates the graph's contours at specified points
  */
  var displayContourFuncArray = function ()
	{
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
    for(var i = 0; i < cArray.length; i++)
    {
      vertexArray[i] = [];
      edgeArray[i] = [];
      
      vertexArray[i] = marchingSquares( xMin, xMax, yMin, yMax, f, cArray[i], 256 );
      
      for (var n = 0; n < vertexArray[i].length; n += 2)
      {
        edgeArray[i].push( vertexArray[i][n], vertexArray[i][n+1] );
      }
      
      edgeData[i].set( "width", edgeArray[i].length/2 );
      edgeData[i].set( "data", edgeArray[i] );
    }
	}
  
	//<!-- CONTOUR CALCULATIONS (MARCHING SQUARES) =============================================== -->
	
	function marchingSquares(xMin, xMax, yMin, yMax, f, c, resolution)
	{
		var xStep = (xMax - xMin) / resolution;
		var yStep = (yMax - yMin) / resolution;
		var points = [];
		for (var x = xMin; x < xMax; x += xStep)
		{
			for (var y = yMin; y < yMax; y += yStep)
			{
				var z1 = f(x,y);				// bottom left corner
				var z2 = f(x+xStep, y);			// bottom right corner
				var z4 = f(x+xStep, y+yStep);	// top right corner
				var z8 = f(x, y+yStep);			// top left corner
				var n = 0;
				if (z1 > c) n += 1;
				if (z2 > c) n += 2;
				if (z4 > c) n += 4;
				if (z8 > c) n += 8;
				
				// calculate linear interpolation values along the given sides.
				//  to simplify, could assume each is 0.5*xStep or 0.5*yStep accordingly.
				var bottomInterp 	= (c - z1) / (z2 - z1) * xStep;
				var topInterp 		= (c - z8) / (z4 - z8) * xStep;
				var leftInterp 		= (c - z1) / (z8 - z1) * yStep;
				var rightInterp 	= (c - z2) / (z4 - z2) * yStep;
				
				// for a visual diagram of cases: https://en.wikipedia.org/wiki/Marching_squares
				if (n == 1 || n == 14) // lower left corner
					points.push( [x, c, y+leftInterp], [x+bottomInterp, c, y] );
					
				else if (n == 2 || n == 13) // lower right corner
					points.push( [x+bottomInterp, c, y], [x+xStep, c, y+rightInterp] );
					
				else if (n == 4 || n == 11) // upper right corner
					points.push( [x+topInterp, c, y+yStep], [x+xStep, c, y+rightInterp] );
					
				else if (n == 8 || n == 7) // upper left corner
					points.push( [x, c, y+leftInterp], [x+topInterp, c, y+yStep] );
					
				else if (n == 3 || n == 12) // horizontal
					points.push( [x, c, y+leftInterp], [x+xStep, c, y+rightInterp] );
					
				else if (n == 6 || n == 9) // vertical
					points.push( [x+bottomInterp, c, y], [x+topInterp, c, y+yStep] );
					
				else if (n == 5) // should do subcase // lower left & upper right
					points.push( [x, c, y+leftInterp], [x+bottomInterp, c, y], [x+topInterp, c, y+yStep], [x+xStep, c, y+rightInterp] );
					
				else if (n == 10) // should do subcase // lower right & upper left
					points.push( [x+bottomInterp, c, y], [x+xStep, c, y+rightInterp], [x, c, y+yStep/2], [x, c, y+leftInterp], [x+topInterp, c, y+yStep] );
					
				else if (n == 0 || n == 15) // no line segments appear in this grid square.
					points.push();
				
			}
		}
		return points;
	
	}
	
	//<!-- END OF CONTOUR CALCULATIONS (MARCHING SQUARES) =============================================== -->
  
  //<!-- START OF CRITICAL PTS CALC  ==================================================== -->
	
	// the function to calculate possible critical points within a set range
	// note this function would return an array of points only
	// it does not result in the points plotted to the graph (yet)
	function critPointsCalc(xMin, xMax, yMin, yMax, f) 
	{	
		 fx = math.derivative(functionText, 'x') ; 
		 fy = math.derivative(functionText, 'y') ; 
		
		var step = 0.01; // density or granularity of search
		var points = [];

		// for every possible x (with increment of "step") within xRange
		for (var x = xMin; x <= xMax; x += step) 
		{	
			// combined with every possible y (with increment of "step") within yRange
			for (var y = yMin; y <= yMax; y += step) 
			{ 
				var z = f(x,y); // calculate the z val of this combination

				let scope = {x:x, y:y}; // to calculate the fx and fy values of this combination
			
				var dzx = fx.evaluate(scope) ;
				var dzdx = Math.round(dzx * 100) / 100; //rounded to remove noise
	
				var dzy = fy.evaluate(scope) ;
				var dzdy = Math.round(dzy * 100) / 100; //rounded to remove noise
				
				var xRound = Math.round(x * 100) / 100;
				var yRound = Math.round(y * 100) / 100;
				var zRound = Math.round(z * 100) / 100;
				
				if (dzdx == 0 && dzdy == 0) {
					points.push( [xRound,zRound,yRound] );
				}
			}
		}
		return points;
	}
	
	//<!-- END OF CRITICAL PTS CALC  ==================================================== -->	
  
 /**
  * updates the graph's critical points at specified points
  */
	var displayCritPointsFunc = function ()
	{
		// create alert first before calculating the critical points, for better UX
		// since calculations will be quite slow, especially for more complex equations
		// further UX improvements can include a loading/buffer screen while calculating
		alert("Critical points will be calculated upon closing this alert.\n\nThis will take a while.\n\nDO NOT refresh the page or press any other buttons while it's loading as it will slow down the process.");
		
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
		var fx = math.derivative(functionText, 'x')  ; 
		var fy = math.derivative(functionText, 'y')  ; 
		
		var fxx = math.derivative(fx, 'x');
		var fyy = math.derivative(fy, 'y');
		var fxy = math.derivative(fx, 'y');
		
		// emptying the array first for every time the function is called
		// so that there's no lingering from previous equation
		critPointsDispArray = [];

		
		// now we transfer the points values from the calculation
		// for this module, we just calculate the critical points within the boundary R
		critPointsDispArray = critPointsCalc( lowerXRange, upperXRange, lowerYRange, upperYRange, f);
		
		critPointData.set( "width", critPointsDispArray.length);
		critPointData.set( "data", critPointsDispArray); // remember that it's x,z,y for plotting the points with Mathbox
		
		// now we are re-indexing the critPointsDispArray
		// for labeling purposes, we want to display the point locations in their proper format x,y (z is not needed)
		critPointsLabelArray = [];
		critPointsClassArray = [];
		
		for (var m = 0; m < critPointsDispArray.length; m++)
		{
			var xPoint = critPointsDispArray[m][0];
			var yPoint = critPointsDispArray[m][2];  // remember that y point is stored in [2] index
			
			critPointsLabelArray.push( [ "(" + xPoint + "," + yPoint + ")" ] );
			
			
			var fxxVal = fxx.evaluate({x: xPoint, y: yPoint}) ;
			var fxx0 = Math.round(fxxVal * 100) / 100;

			var fyyVal = fyy.evaluate({x: xPoint, y: yPoint}) ;
			var fyy0 = Math.round(fyyVal * 100) / 100;

			var fxyVal = fxy.evaluate({x: xPoint, y: yPoint}) ;
			var fxy0 = Math.round(fxyVal * 100) / 100;

			var H = fxx0 * fyy0 - Math.pow(fxy0,2);

			if (H == 0) {
				critPointsClassArray.push( "undefined" );
			} else if (H < 0) {
				critPointsClassArray.push( "saddle point" );
				
			} else if (H > 0) {
				if (fxx0 < 0) {
					critPointsClassArray.push( "maxima" );

				} else if (fxx0 > 0) {
					critPointsClassArray.push( "minima" );

				}
			}	
		}
		
		critPointLocText.set( "width", critPointsLabelArray.length);

		critPointClassText.set( "width", critPointsClassArray.length);

		/*globalMaximaArray = [];
		globalMaximaData.set("data",[globalMaximaArray]);
		
		globalMinimaArray = [];
		globalMinimaData.set("data",[globalMinimaArray]);
		
		globalMaximaLocArray = [];
		globalMinimaLocArray = [];
		
		
		optimaArray = [];
		optimaArray = OptimaCalcFunc(lowerXRange, upperXRange, lowerYRange, upperYRange, f);*/
		
		//OptimaDisplayFunc();

	}
  
 /**
  * updates the graph's critical points at specified points
  */
  var displayTraceFunc = function ()
	{
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );

		traceZ = f(traceX, traceY); // z0 = f(x0,y0)
		tracePointData.set("data", [ [traceX, traceZ, traceY] ] );
		
		var traceXDisp = Math.round(traceX * 100) / 100;
		var traceYDisp = Math.round(traceY * 100) / 100;
		var traceZDisp = Math.round(traceZ * 100) / 100;
		traceText.set( "data", ["(" + traceXDisp + " , " + traceYDisp + " , " + traceZDisp + ")"] );
		traceTextView.set("visible", traceVisible);
		

	}
  
 /**
  * updates the range for searching global optima
  */  
	var displayRangeFunc = function()
	{
		f = Parser.parse( functionText ).toJSFunction( ['x','y'] );
		
		xMinCurveData.set("expr", 
		    function (emit, y, j, t) 
			{ 
			    emit( lowerXRange, f(lowerXRange,y), y );
			}
		);
		xMinCurveText.set("data", ["xMin = " + Math.round(lowerXRange * 100) / 100]);

		xMaxCurveData.set("expr", 
		    function (emit, y, j, t) 
			{ 
			    emit( upperXRange, f(upperXRange,y), y );
			}
		);
		xMaxCurveText.set("data", ["xMax = " + Math.round(upperXRange * 100) / 100]);
		
		yMinCurveData.set("expr", 
		    function (emit, x, i, t) 
			{ 
			    emit( x, f(x,lowerYRange), lowerYRange );
			}
		);
		yMinCurveText.set("data", ["yMin = " + Math.round(lowerYRange * 100) / 100]);
		
		yMaxCurveData.set("expr", 
		    function (emit, x, i, t) 
			{ 
			    emit( x, f(x,upperYRange), upperYRange );
			}
		);
		yMaxCurveText.set("data", ["yMax = " + Math.round(upperYRange * 100) / 100]);
	}
  
  //OPTIMA DISPLAY FUNCTIONS-----------------------------------------------------------------
  // first we find possible z values of every x & y within the Range
	var OptimaCalcFunc = function(xMin, xMax, yMin, yMax, f)
	{
		var points = [];
		var step = 1;
	
		for (var x=xMin; x <= xMax; x += step)
		{
			for (var y=yMin; y <= yMax; y += step)
			{
				var z = f(x,y);
				points.push([x,z,y]);
			}
		}
		return points;
	}
  
  
	var OptimaDisplayFunc = function()
	{
		
		//we run through the z values in optimaArray first
		var zPointsArray = [];
		for (var u = 0; u < optimaArray.length; u++) 
		{
			zPointsArray.push(optimaArray[u][1]);
		}
		
		//then we find what is the biggest Z value and the smallest Z value
		var max = Math.max(...zPointsArray);
		var min = Math.min(...zPointsArray);
			
		for (var v = 0; v < optimaArray.length; v++) 
		{	
			var globalX = Math.round(optimaArray[v][0] * 100) / 100;
			var globalY = Math.round(optimaArray[v][2] * 100) / 100;
			
			if (optimaArray[v][1] == max)
			{
				globalMaximaArray.push(optimaArray[v]);
				globalMaximaLocArray.push([ "(" + globalX + "," + globalY + ")" ]);
			}
			if (optimaArray[v][1] == min)
			{
				globalMinimaArray.push(optimaArray[v]);
				globalMinimaLocArray.push([ "(" + globalX + "," + globalY + ")" ]);
			}
		}
		
		globalMaximaData.set("width",globalMaximaArray.length);
		globalMaximaData.set("data",globalMaximaArray);
		
		globalMaximaText.set("width",globalMaximaArray.length);
		globalMaximaLocText.set("width",globalMaximaArray.length);
		
		
		globalMinimaData.set("width",globalMinimaArray.length);
		globalMinimaData.set("data",globalMinimaArray);
		
		globalMinimaText.set("width",globalMinimaArray.length);
		globalMinimaLocText.set("width",globalMinimaArray.length);
		
	}
	//OPTIMA DISPLAY FUNCTIONS END-------------------------------------------------------------