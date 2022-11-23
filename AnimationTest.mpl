AnimationTest := module()
    option package;

    local diffEqs:

    export NumberGenerator:
    
    export mlaPath := cat(kernelopts(homedir),"/maple/Packages/AnimationTest/AnimationTest.mla"):

    module NumberGenerator()
	    option object;

	    local ModuleApply, uniformStart, uniformEnd, nSamples, samples, sampleTimes, seed, rNumbers, currentSamplesIndex;
		
		# Called when NumberGenerator is called procedure-style
		ModuleApply::static := proc()
			Object(NumberGenerator, _passed);
		end;

		# Called when NumberGenerator is instantiated using Object
		export ModuleCopy::static := proc( self::NumberGenerator, proto::NumberGenerator, uniformStart, uniformEnd, nSamples, seed, $)
			# If no extra arguments other then self and proto
			if ( _npassed = 2 ) then
				self:-uniformStart := proto:-uniformStart;
				self:-uniformEnd := proto:-uniformEnd;
				self:-nSamples := proto:-nSamples;
				self:-seed := proto:-seed;
				self:-samples := generateSamples(uniformStart, uniformEnd, nSamples,seed);
				self:-currentSamplesIndex := proto:-currentSamplesIndex;
				self:-sampleTimes := proto:-sampleTimes;
				self:-rNumbers := proto:-rNumbers;
			else
				self:-uniformStart := uniformStart;
				self:-uniformEnd := uniformEnd;
				self:-nSamples := nSamples;
				self:-seed := seed;
				self:-samples := generateSamples(uniformStart, uniformEnd, nSamples,seed):
				self:-currentSamplesIndex := 1;
				self:-sampleTimes := Array([]);
				self:-rNumbers := Table([]);
			end;
		end;

	    export generate := proc(t)
			local currentSize:
			currentSize := ArrayTools:-Size(sampleTimes)[2]:
			sampleTimes(currentSize+1) := t:
	    	if (not assigned(rNumbers[t])) then
				rNumbers[t] := samples[currentSamplesIndex];
				currentSamplesIndex := currentSamplesIndex + 1;
			end; 

			return rNumbers[t];
		end;
		
		export generateSamples::static := proc(uniformStart, uniformEnd, nSamples,seed, $)
			setState(seed);
			return Statistics:-Sample('Uniform'(uniformStart, uniformEnd))(nSamples);
		end;
		
		export reset := proc(uniformStart,uniformEnd, nSamples, newSeed) 
			samples := generateSamples(uniformStart, uniformEnd, nSamples,seed):
			seed := seed;
			currentSamplesIndex := 1; 
		end;
		
		export viewSamples::static := proc(self::NumberGenerator) 
			print(self:-samples); 
		end;

		export viewTable::static := proc(self::NumberGenerator)
			print(self:-rNumbers);
		end;

		export viewGenerated::static := proc(self::NumberGenerator)
			local i, currentSize:
			currentSize := ArrayTools:-Size(self:-sampleTimes)[2]:
			[seq([self:-sampleTimes[i], self:-samples[i]], i = 1..currentSize)]:
		end:

		export getAtTime := proc(t)
		  return rNumbers[t];
		end;
	end:

    diffEqs := proc()
        diff(I*f(t),t) = (1+'gen(t)')*f(t)/2+(0.1*cos(t))*g(t),
        diff(I*g(t),t) = -(1+'gen(t)')*g(t)/2+(0.1*cos(t))*f(t);
    end:

    export solveEqs := proc(seedn, step, Np)
        local p, i, Sz, Sx, Sy, slist, eqsyscomplex, fcns, ng:
        global gen;

        ng := NumberGenerator(-1,1, 10000000, seedn);

        gen := proc(x) 
            if type(x, numeric) then 
                return generate(x,ng):
            end:

		    return 'gen(t)':
	    end;

        fcns := {f(t),g(t)}:

        eqsyscomplex := diffEqs():

        print(eqsyscomplex, seedn, step, Np, test());

        p := dsolve(
            {eqsyscomplex, f(0)=1, g(0)=0}, 
            fcns, 
            type=numeric, 
            method=rkf45, 
            maxfun=0, 
            abserr=0.0001, 
            relerr=0.0001, 
            output=Array([seq(i/step, i = 0 .. Np)])):

        Sz := seq(abs(p(2)[i,2])^2-abs(p(2)[i,3])^2, i=1..Np):
		Sx := seq(conjugate(p(2)[i,2])*p(2)[i,3]+conjugate(p(2)[i,3])*p(2)[i,2], i=1..Np):
		Sy := seq(-I*conjugate(p(2)[i,2])*p(2)[i,3]+I*conjugate(p(2)[i,3])*p(2)[i,2],i=1..Np):

		slist := [seq([-Re(Sx[i]), -Re(Sy[i]), -Re(Sz[i])], i = 1 .. Np)]:

        return slist;
    end;

    export computeAverage := proc(step, Np, N)
		local l, i, S;   


        S := Grid:-Map(solveEqs, [seq(i, i = 1 .. N)], step, Np):      
		
		l := S[1]:    # first list in S          
		
		for i from 2 to N do           
			l:= l+S[i];     # iteratively sum the lists in S element-wise        
		end do:          
		
		l:= l/N;     # divide the elements element-wise to final list by N to get average      

		return l;
	end;

	export animate := proc(step, Np, N)
		local Stot, q1, animationData;

		Stot := computeAverage(step, Np, N):

		q1 := Matrix(Stot):

		animationData := computeAnimationData(q1):

		displayAnimation(seq(animationData));
	end;

	export computeAnimationData := proc(q1)
		local windowSize, parrow, ptrail, finalTime, times, psphere, m1, m2;

		m1 := proc(q1, t, windowSize) 
			# for the number of points less than or equal to windowSize, the window size equals the number of points
			if t <= windowSize then 
				return plots:-pointplot3d(q1[1 .. t, () .. ()], connect = true); 
			end if; 
			
			# after the windowSize-th point, the window size is fixed at windowSize
			return plots:-pointplot3d(q1[t - windowSize + 1 .. t, () .. ()], connect = true); 
		end;

		m2 := proc(q1, t) return plots:-arrow([[0, 0, 0], q1[t]]); end;

		finalTime := LinearAlgebra:-Dimensions(q1)[1];
		times := [seq(1 .. finalTime, 1)];
		windowSize := 1000;
		parrow := plots:-animate(m2, [q1, t], t = times, view = [-1 .. 1, -1 .. 1, -1 .. 1], color = red);
		ptrail := plots:-animate(m1, [q1, t, windowSize], t = [seq(1 .. finalTime)], color = red);
		psphere := plottools:-sphere([0, 0, 0], 1, transparency = 0.9);
		return [ptrail, parrow, psphere]:	
	end:	

    export displayAnimation := proc(ptrail, parrow, psphere)
		plots:-display([ptrail, parrow, psphere], scaling = constrained, axes = boxed);
	end: 
end:

LibraryTools:-Save('AnimationTest', mlaPath);