Nneurons = [10,10]; % the number of neurons in each sequence
Dt = [10,10]; % the number of time steps between each neuron in the sequence
Pseq = [1,1]; % the probability of the sequence occuring
NeuronNoise = 0.001; % the noise in a neurons firing rate
SeqNoiseTime = [0,0]; % the noise in the sequence aka jitter (p of each neuron jittered 1 dt)
SeqNoiseNeuron = [1,1]; % the probability that a neuron participates in a given seq
T = 10000;
shared  = 1;
%%
V = data.sequences(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,shared);