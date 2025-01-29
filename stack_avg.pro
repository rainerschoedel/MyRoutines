; $Id: stack_median.pro, v 1.0 Aug 1999 e.d. $
;
;+
; NAME:
;	STACK_AVG
;
; PURPOSE:
;	Combine frames in a 3D stack by median superposition.
;       This is a modified version of StarFinder STACK_MEDIAN
;       that lets the user choose between median and mean of a stack.
; CATEGORY:
;	Signal processing.
;
; CALLING SEQUENCE:
;	Result = STACK_MEDIAN(Stack)
;
; INPUTS:
;	Stack:	3D array, representing the stack of frames to be combined
;
; KEYWORD PARAMETERS:
;	MASK:	3D binary array used to mask the frames in the Stack.
;		The n-th frame of the 3D array Mask passed with this keyword may
;		be defined as follows:
;		Mask[j, i, n] = 1, if  Stack[j, i, n] is a valid pixel
;		              = 0, if  Stack[j, i, n] must be rejected
;		The default is no masking.
;
;	WEIGHTS:	Set this keyword to a N-components vector, where N is the
;		number of frames in the Stack. The n-th component is the weight to
;		apply to the n-th input frame.
;		In practice the weights are divided by their minimum value and
;		rounded to integers. Every plane in the Stack is then replicated a
;		number of times equal to the corresponding weight.
;
; 
;       AVGYN: If AVGYN = 0 (default), then the median of a stack is computed, else the average.
;
;       CLIP: If set, then the CLIP highest and lowest values in the stack
;             will be ignored.
; 
; OUTPUTS:
;	Result:	2D array, representing the median of the Stack planes
;
; MODIFICATION HISTORY:
; 	Written by:	Emiliano Diolaiti, August 1999.
;       Modified by:    Rainer Schoedel, January 2012, September 2012,
;                       January 2018

FUNCTION stack_avg, stack, MASK = mask, WEIGHTS = weights, AVGYN = AVGYN

        if not(keyword_set(AVGYN)) then AVGYN = 0
	on_error, 2
	if  size52(stack, /N_DIM) ne 3  then  return, stack
	s = size52(stack, /DIM)  &  s0 = s[0]  &  s1 = s[1]  &  n_frames = s[2]
	; Mask some frames?
	masked = n_elements(mask) ne 0
	; Weighted median?
	weighted = n_elements(weights) ne 0
	if  weighted  then  weighted = min(weights) gt 0
	if  weighted  then begin
	   w = weights > 0  &  w = round(w/min(w))
	endif else  w = replicate(1, n_frames)
	; Compute median
	m = make_array(s0, s1, TYPE = size52(stack, /TYPE))
	for  i = 0L, s1 - 1  do  for  j = 0L, s0 - 1  do begin
	   slice = stack[j,i,*]  &  w_ji = w  &  n = n_frames
	   ; Mask undesired frames
	   if  masked  then begin
	      accept = where(mask[j,i,*] ne 0, n)
	      if  n ne 0  then begin
	         slice = slice[accept]  &  w_ji = w_ji[accept]
	      endif
	   endif
	   ; Replicate frames according to their weight
	   if  weighted and n ne 0  then begin
	      temp = replicate(slice[0], w_ji[0])
	      for  k = 1L, n - 1  do $
	         temp = append_elements(temp, replicate(slice[k], w_ji[k]))
	      slice = temp
	   endif
	   ; Compute median of pixel [j,i]
	   if  (n ne 0)  then if (avgyn ne 0) then  begin
            RESISTANT_Mean,slice,3.0,mean_slice,sigma_bg,Num_Rej
            m[j,i] = mean_slice
          endif else begin
            m[j,i] = median(slice, /EVEN)
          endelse
	endfor
      return, m
end
