function check_aliasing(d_frep, T_batch, NT)
    d_frep = 10e6;        % 10 MHz
    comb2_shift = 0;      % Let's set to 0 for a moment to avoid Trap 1
    T_batch = 4e-9;       % 4 ns batch
    NT = 100;             % Look at 100 batches
    
    % 2. Pick a comb line
    n = 1;
    delta_n = (n * d_frep) - comb2_shift; 
    
    % 3. Calculate the phase step per batch
    cycles_per_batch = delta_n * T_batch;
    disp(['Cycles per batch: ', num2str(cycles_per_batch)]);
    
    % 4. Generate the phasors
    j = 1:NT;
    Tj = (j-1) * T_batch;
    phasors = exp(1i * 2 * pi * delta_n * Tj);
    
    % 5. Plot the circle!
    figure;
    plot(real(phasors), imag(phasors), 'o-', 'LineWidth', 1.5);
    grid on; axis equal;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    title(sprintf('Phasor Rotation for n=%d', n));
    xlabel('Real'); ylabel('Imaginary');
end