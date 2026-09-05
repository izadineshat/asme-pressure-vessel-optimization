% =========================================================================
% PROJECT:     TECHNO-ECONOMIC ASME SEC VIII DIV 1 VESSEL OPTIMIZATION
% Author:        Reza Neshat
% Affiliation:    Pars Mabna Energy | Mechanical Engineering (Applied Design)
% Contact:       izadi@pars-mabna.com | izadineshat@gmail.com
% BENCHMARK:   SANDGREN (1990) / KANNAN & KRAMER (1994) / ASME BPVC UG-27/32
% CURRENCY:    IRANIAN RIAL (IRR)
% SCRIPT:      run_iran_asme_vessel_optimization.m
% =========================================================================

function [best_sol, best_cost_irr, conv_history] = run_iran_asme_vessel_optimization()
    clc;
    close all;
    
    % Deterministic Pseudo-Random Seed
    rng(100, 'twister');
    
    % Commercial plate gauge catalog available in Iran (DIN EN 10029 / A516-Gr70)
    plate_catalog = [6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60]; % mm
    
    % Algorithm Hyperparameters (Constrained Vectorized PSO)
    swarmsize = 120;
    max_iter  = 400;
    w_start   = 0.90;
    w_end     = 0.35;
    c1        = 1.75;
    c2        = 1.95;
    
    % Design Bounds: [Ts (mm), Th (mm), R (mm), L (mm)]
    lb = [min(plate_catalog), min(plate_catalog),  500.0,  1000.0];
    ub = [max(plate_catalog), max(plate_catalog), 3000.0,  8000.0];
    dim = 4;
    
    v_max = 0.15 * (ub - lb);
    v_min = -v_max;
    
    % Swarm Initialization
    X = repmat(lb, swarmsize, 1) + repmat((ub - lb), swarmsize, 1) .* rand(swarmsize, dim);
    V = repmat(v_min, swarmsize, 1) + repmat((v_max - v_min), swarmsize, 1) .* rand(swarmsize, dim);
    
    % Discrete Snapping to Commercial Plate Catalog
    X(:, 1) = snap_to_catalog(X(:, 1), plate_catalog);
    X(:, 2) = snap_to_catalog(X(:, 2), plate_catalog);
    
    P_best = X;
    P_best_fitness = zeros(swarmsize, 1);
    
    for i = 1:swarmsize
        [fit, ~, ~] = evaluate_fitness(X(i, :), 1);
        P_best_fitness(i) = fit;
    end
    
    [G_best_fitness, g_idx] = min(P_best_fitness);
    G_best = P_best(g_idx, :);
    
    conv_history = zeros(max_iter, 1);
    raw_cost_history = zeros(max_iter, 1);
    
    % Iterative Optimization Loop
    for t = 1:max_iter
        % Non-linear dynamic inertia reduction
        w = w_end + (w_start - w_end) * ((max_iter - t) / max_iter)^1.2;
        
        r1 = rand(swarmsize, dim);
        r2 = rand(swarmsize, dim);
        
        V = w * V + c1 * r1 .* (P_best - X) + c2 * r2 .* (repmat(G_best, swarmsize, 1) - X);
        V = max(min(V, repmat(v_max, swarmsize, 1)), repmat(v_min, swarmsize, 1));
        
        X = X + V;
        
        % Box Boundary Handling
        X = max(min(X, repmat(ub, swarmsize, 1)), repmat(lb, swarmsize, 1));
        
        % Continuous-to-Discrete Projection on Physical Steel Gauges
        X(:, 1) = snap_to_catalog(X(:, 1), plate_catalog);
        X(:, 2) = snap_to_catalog(X(:, 2), plate_catalog);
        
        for i = 1:swarmsize
            [fit_val, ~, ~] = evaluate_fitness(X(i, :), t);
            
            if fit_val < P_best_fitness(i)
                P_best(i, :) = X(i, :);
                P_best_fitness(i) = fit_val;
            end
        end
        
        [current_best_fit, min_idx] = min(P_best_fitness);
        if current_best_fit < G_best_fitness
            G_best_fitness = current_best_fit;
            G_best = P_best(min_idx, :);
        end
        
        conv_history(t) = G_best_fitness;
        raw_cost_history(t) = compute_base_cost(G_best);
        
        if mod(t, 50) == 0 || t == max_iter
            fprintf('Iter: %03d | Penalized Cost: %16.0f IRR | True Cost: %16.0f IRR\n', ...
                t, G_best_fitness, raw_cost_history(t));
        end
    end
    
    best_sol = G_best;
    best_cost_irr = compute_base_cost(best_sol);
    
    % Terminal Audit Output
    print_industrial_audit(best_sol, best_cost_irr);
    
    % Diagnostics Generation & High-Resolution Export
    generate_and_export_diagnostics(conv_history, raw_cost_history, best_sol);
end

function [fitness, raw_cost, max_viol] = evaluate_fitness(x, iter)
    raw_cost = compute_base_cost(x);
    g = compute_asme_constraints(x);
    
    violations = max(0, g);
    max_viol = max(violations);
    
    % Dynamic Exterior Quadratic Penalty Formulation
    gamma = 1.2e8;
    beta  = 1.6;
    penalty_factor = (gamma * (iter / 100))^beta;
    
    % Non-dimensional constraint normalization vector
    norm_viol = [violations(1) / 50.0;
                 violations(2) / 50.0;
                 violations(3) / 2.5e10;
                 violations(4) / 8000.0];
                 
    penalty = penalty_factor * sum(norm_viol.^2);
    fitness = raw_cost + penalty;
end

function cost_irr = compute_base_cost(x)
    Ts = x(1);
    Th = x(2);
    R  = x(3);
    L  = x(4);
    
    % Cost calibration factors for Iranian Industrial Sector (IRR)
    C1 = 36.98;   % Rolled shell material (750,000 IRR/kg)
    C2 = 118.32;  % Head material and dishing operations (1,200,000 IRR/kg)
    C3 = 35.00;   % Longitudinal SAW seam weld
    C4 = 263.89;  % Circumferential girth joint weld
    
    cost_irr = C1 * Ts * R * L + ...
               C2 * Th * (R^2) + ...
               C3 * (Ts^2) * L + ...
               C4 * (Ts^2) * R;
end

function g = compute_asme_constraints(x)
    Ts = x(1);
    Th = x(2);
    R  = x(3);
    L  = x(4);
    
    g = zeros(4, 1);
    % g1: Shell Hoop Stress Requirement (UG-27)
    g(1) = (0.02994 * R + 2.0) - Ts;
    
    % g2: Head Membrane Stress Requirement (UG-32)
    g(2) = (0.01475 * R + 2.0) - Th;
    
    % g3: Critical Internal Capacity Volume (V >= 25 m^3)
    V_internal = pi * (R^2) * L + (4.0 / 3.0) * pi * (R^3);
    g(3) = 2.5e10 - V_internal;
    
    % g4: Transportation Envelope Constraint
    g(4) = L - 8000.0;
end

function snapped_vals = snap_to_catalog(vals, catalog)
    snapped_vals = zeros(size(vals));
    for i = 1:length(vals)
        [~, min_id] = min(abs(catalog - vals(i)));
        snapped_vals(i) = catalog(min_id);
    end
end

function print_industrial_audit(sol, cost_irr)
    g = compute_asme_constraints(sol);
    
    Ts = sol(1); Th = sol(2); R = sol(3); L = sol(4);
    C1 = 36.98; C2 = 118.32; C3 = 35.00; C4 = 263.89;
    
    cost_shell_mat = C1 * Ts * R * L;
    cost_head_mat  = C2 * Th * (R^2);
    cost_long_weld = C3 * (Ts^2) * L;
    cost_circ_weld = C4 * (Ts^2) * R;
    
    v_net = (pi * (R^2) * L + (4.0/3.0)*pi*(R^3)) * 1e-9;
    
    fprintf('\n=======================================================================\n');
    fprintf('  PARS MABNA ENERGY - ASME SEC VIII DIV 1 OPTIMIZATION AUDIT SHEET     \n');
    fprintf('=======================================================================\n');
    fprintf('Optimal Geometrical Parameters:\n');
    fprintf('  Shell Plate Thickness (Ts): %6.1f mm [Standard Commercial Sheet]\n', Ts);
    fprintf('  Head Plate Thickness (Th) : %6.1f mm [Standard Commercial Sheet]\n', Th);
    fprintf('  Internal Radius (R)       : %8.2f mm (ID: %.2f mm)\n', R, 2*R);
    fprintf('  Tangent-to-Tangent (L)    : %8.2f mm\n', L);
    fprintf('  Total Volume Capacity     : %8.3f m^3\n', v_net);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('Financial Breakdown (Iranian Rial):\n');
    fprintf('  Shell Material Purchase   : %16.0f IRR (%.1f%%)\n', cost_shell_mat, (cost_shell_mat/cost_irr)*100);
    fprintf('  Head Material & Forming   : %16.0f IRR (%.1f%%)\n', cost_head_mat, (cost_head_mat/cost_irr)*100);
    fprintf('  Longitudinal Seam Welding : %16.0f IRR (%.1f%%)\n', cost_long_weld, (cost_long_weld/cost_irr)*100);
    fprintf('  Circumferential Girth Weld: %16.0f IRR (%.1f%%)\n', cost_circ_weld, (cost_circ_weld/cost_irr)*100);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('TOTAL FABRICATED COST       : %16.0f IRR\n', cost_irr);
    fprintf('TOTAL EQUIVALENT (TOMAN)    : %16.0f Tomans\n', cost_irr / 10);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('ASME Constraint Margins [g_i(x) <= 0]:\n');
    fprintf('  g1 (Shell Hoop Stress)    : %12.4f mm   | Status: %s\n', g(1), check_feasibility(g(1)));
    fprintf('  g2 (Head Membrane Stress) : %12.4f mm   | Status: %s\n', g(2), check_feasibility(g(2)));
    fprintf('  g3 (Deficit Volume)       : %12.4e mm^3 | Status: %s\n', g(3), check_feasibility(g(3)));
    fprintf('  g4 (Length Exceedance)    : %12.4f mm   | Status: %s\n', g(4), check_feasibility(g(4)));
    fprintf('=======================================================================\n\n');
end

function str = check_feasibility(val)
    if val <= 1e-4
        str = 'SATISFIED (COMPLIANT)';
    else
        str = 'VIOLATED';
    end
end

function generate_and_export_diagnostics(conv_curve, raw_cost_curve, sol)
    fig = figure('Color', [1 1 1], 'Position', [100, 100, 1200, 800]);
    
    Ts = sol(1); Th = sol(2); R = sol(3); L = sol(4);
    
    % Panel 1: Convergence History
    subplot(2, 2, 1);
    plot(1:length(conv_curve), conv_curve / 1e9, 'r-', 'LineWidth', 1.2); hold on;
    plot(1:length(raw_cost_curve), raw_cost_curve / 1e9, 'b-', 'LineWidth', 1.8);
    grid on; box on;
    xlabel('Iteration Index', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Cost (Billion IRR)', 'FontSize', 10, 'FontWeight', 'bold');
    title('Convergence History: Penalized vs True Structural Cost', 'FontSize', 11);
    legend('Penalized Fitness', 'Real Structural Cost', 'Location', 'northeast');
    
    % Panel 2: Cost Breakdown
    subplot(2, 2, 2);
    C1 = 36.98; C2 = 118.32; C3 = 35.00; C4 = 263.89;
    costs = [C1*Ts*R*L, C2*Th*(R^2), C3*(Ts^2)*L, C4*(Ts^2)*R];
    labels = {'Shell Plate', 'Head Plate & Forming', 'Longitudinal Weld', 'Circumferential Weld'};
    p = pie(costs, labels);
    title('Fabrication Cost Distribution (IRR)', 'FontSize', 11);
    
    % Panel 3: Active Boundary Mapping
    subplot(2, 2, 3);
    r_sweep = linspace(500, 2500, 250);
    t_req = 0.02994 * r_sweep + 2.0;
    plot(r_sweep, t_req, 'k--', 'LineWidth', 1.5); hold on;
    plot(R, Ts, 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
    grid on; box on;
    xlabel('Internal Radius R (mm)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Shell Thickness Ts (mm)', 'FontSize', 10, 'FontWeight', 'bold');
    title('Design Feasibility vs ASME UG-27 Minimum Boundary', 'FontSize', 11);
    legend('ASME UG-27 Threshold', 'Optimal Commercial Point', 'Location', 'northwest');
    
    % Panel 4: Constraint Residual / Slack Analysis
    subplot(2, 2, 4);
    g = compute_asme_constraints(sol);
    bar_vals = [g(1), g(2), g(3)/1e9, g(4)/100];
    b = bar(bar_vals, 'FaceColor', [0.15 0.45 0.68]);
    grid on; box on;
    set(gca, 'XTickLabel', {'g_1 (Shell)', 'g_2 (Head)', 'g_3 (Vol \times 10^{-9})', 'g_4 (L / 100)'});
    ylabel('Constraint Slack Value', 'FontSize', 10, 'FontWeight', 'bold');
    title('ASME Constraint Active Status (Slack Analysis)', 'FontSize', 11);
    
    % High-Resolution Print Export (300 DPI)
    output_filename = 'asme_vessel_optimization_diagnostics.png';
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, output_filename, '-dpng', '-r300');
    fprintf('Diagnostic figure successfully exported at 300 DPI: %s\n', output_filename);
end