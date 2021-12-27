function stop = myOutputFcn(optimValues,state)
        stop = false;
        switch state
%             case 'iter'
                case 'done'
                
                % Make updates to plot or guis as needed
                
                %% PRINT RESULTS
                fprintf(' It.:%5i Obj.:%11.4f ch.:%7.3f\n',optimValues.iteration,optimValues.fval,optimValues.stepsize);
                %% PLOT DENSITIES
                
            otherwise
        end % switch
    end % myOutputFcn