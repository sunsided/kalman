# Unscented Kalman Filter for Frequency Tracking

This example uses a Scaled Unscented Kalman Filter in order to track a signal of the form

	y(t) = a*sin(b*t+c) + cos(d*t) + e

using a model of the form

	z(t) = a*sin(b*t+c) + d

that is, the state vector consists of frequency, phase, amplitude and offset. 