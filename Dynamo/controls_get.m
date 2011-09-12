function ret = controls_get(control_mask)
global OC;

if nargin == 0
    ret = OC.seq.raw_controls;
else
    ret = OC.seq.raw_controls(control_mask);
end

