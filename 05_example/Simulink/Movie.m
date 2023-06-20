classdef Movie < handle
	%MOVIE_OPENSOURCE This class creates a user interface than enables to
	%visually analyse the simulated output trajectories
	
	properties
		frame_idx	= 1;
		play		= true;
		F
		f
		ax
		b
		s
	end
	
	methods
		
		function obj = Movie(F)
			%MOVIE_OPENSOURCE Construct an instance of this class

			% Movie frames
			obj.F	= F;
			% Figure
			obj.f	= uifigure(...
				'Position',[284 257 869 500],...
				'Color','w',...
				'MenuBar','none',...
				'NumberTitle','off',...
				'Name','Simulation Outputs');
			set(obj.f,'CloseRequestFcn',@(src,event) onClose(obj,src,event));
			% Axis
			obj.ax	= uiaxes(obj.f,...
				'Position',[-35 -20 955 560],...
                'BackGroundColor','w');
			% Button
			obj.b	= uibutton(obj.f,...
				'Position',[200 15 150 25],...
				'Text','Pause',...
				'ButtonPushedFcn',@(handle,event) PlayButtonPushed(obj,handle,event));
			% Slider
			obj.s	= uislider(obj.f,...
				'Position',[450 40 300 3],...
				'Limits',[1 numel(F)],...
				'MajorTicks',round(linspace(1,numel(F),numel(F)/2/5+1)),...
				'MajorTickLabels',string(0 : 5 : floor(numel(F))/2),...
				'Value',1,...
				'ValueChangingFcn',@(handle,event) SliderChanging(obj,handle,event),...
				'ValueChangedFcn',@(handle,event) SliderChanged(obj,handle,event));			
		end
		
		%% PlayButton
		function obj = PlayButtonPushed(obj,handle,event)
            obj.play = ~obj.play;
			switch obj.b.Text
				case 'Pause'
					obj.b.Text = 'Play';
				case 'Play'
					obj.b.Text = 'Pause';
			end
		end
		
		%% SliderChanging
		function obj = SliderChanging(obj,handle,event)
			obj.play	= false;
			obj.b.Text	= 'Play';
			obj.frame_idx	= round(event.Value);
		end

		%% SliderChanged
		function obj = SliderChanged(obj,handle,event)
			obj.frame_idx = round(handle.Value);
		end
		
		%% PlayMovie
		function obj = PlayMovie(obj)
            % Ensure movie starts at beginning
            drawnow;
            % Visualize trajectories
			while isvalid(obj)
                if ishandle(obj.ax)
                    imshow(obj.F(obj.frame_idx).cdata,'Parent',obj.ax);
                    obj.s.Value = obj.frame_idx;
                end
				if obj.play
					if obj.frame_idx == numel(obj.F)
						obj.frame_idx = 1;
					else
						obj.frame_idx = obj.frame_idx+1;
					end
				end
				pause(0.1);
			end
		end
		
		%% Close Figure
        function onClose(obj,src,event)
            delete(src)
            delete(obj)
        end
		
	end
end
