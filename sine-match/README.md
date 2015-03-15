# Unscented Kalman Filter for Frequency Tracking

This example uses a additive-noise scaled unscented Kalman filter in order to track a signal of the form

	y(t) = a*sin(b + 2*pi*c*T)

that is, the state vector consists of frequency, phase, amplitude and offset. 