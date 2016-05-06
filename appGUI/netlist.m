function [s] = netlist(x, filename)
    outfile = fopen(filename, 'w');
    fprintf(outfile,'*** MODEL DESCRIPTION ***\n');
    fprintf(outfile,'.model mod1 npn bf=150 \n');
    fprintf(outfile,'.model mod2 pnp bf=150 \n');
    fprintf(outfile,'\n*** NETLIST DESCRIPTION ***\n');
    counts = zeros(1,6); %1 = source, 2 = cap, 3 = ind, 4 = mos, 5 = pnp, 6 = resistor
    for ii = [1:size(x,1)]
        switch x(ii,1)
            case 1  %AC source, value field = voltage, node2 = frequency
                counts(1,1) = counts(1,1)+1;
                fprintf(outfile,'V')
                fprintf(outfile,'%d',counts(1,1));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3))
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4))
                fprintf(outfile,' ')
                fprintf(outfile,'sin(0 %.2f %.2f 0 0)',x(ii,2),x(ii,5));
            case 2 %Capacitor
                counts(1,2) = counts(1,2)+1;
                fprintf(outfile,'C')
                fprintf(outfile,'%d',counts(1,2));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3))
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4))
                fprintf(outfile,' ')
                fprintf(outfile,'%.1e',x(ii,2))
            case 3 %DC source
                counts(1,1) = counts(1,1)+1;
                fprintf(outfile,'V')
                fprintf(outfile,'%d',counts(1,1));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3))
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4))
                fprintf(outfile,' ')
                fprintf(outfile,'dc')
                fprintf(outfile,' ')
                fprintf(outfile,'%.2f',x(ii,2))
            case 4 %Ground, ignore this?
%                 fprintf(outfile,'GND')
%                 fprintf(outfile,'%d',counts(1,component));
            case 5 %Inductor
                counts(1,3) = counts(1,3)+1;
                fprintf(outfile,'L')
                fprintf(outfile,'%d',counts(1,3));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3))
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4))
                fprintf(outfile,' ')
                fprintf(outfile,'%.1f',x(ii,2))
            case 6 %NMOS, node0 = drain, node1 = gate, node2 = source
                counts(1,4) = counts(1,4)+1;
                fprintf(outfile,'M')
                fprintf(outfile,'%d',counts(1,4));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,5));
                fprintf(outfile,' ')
                fprintf(outfile,'mod3');
            case 7 %NPN, node0 = collector, node1 = base, node2 = emitter
                counts(1,5) = counts(1,5)+1;
                fprintf(outfile,'Q')
                fprintf(outfile,'%d',counts(1,5));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,5));
                fprintf(outfile,' ')
                fprintf(outfile,'mod1');
            case 8 %PMOS, node0 = drain, node1 = gate, node2 = source
                counts(1,4) = counts(1,4)+1;
                fprintf(outfile,'M')
                fprintf(outfile,'%d',counts(1,4));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,5));
                fprintf(outfile,' ')
                fprintf(outfile,'mod4');
            case 9 %PNP, node0 = collector, node1 = base, node2 = emitter
                counts(1,5) = counts(1,5)+1;
                fprintf(outfile,'Q')
                fprintf(outfile,'%d',counts(1,5));
                fprintf(outfile,'%d',x(ii,3));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,5));
                fprintf(outfile,' ')
                fprintf(outfile,'mod2');
            case 10 %Resistor
                counts(1,6) = counts(1,6)+1;
                fprintf(outfile,'R')
                fprintf(outfile,'%d',counts(1,6));
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,3))
                fprintf(outfile,' ')
                fprintf(outfile,'%d',x(ii,4))
                fprintf(outfile,' ')
                fprintf(outfile,'%.1f',x(ii,2))
            otherwise
                fprintf(outfile,'Error')
        end
        fprintf(outfile,'\n')
    end
    fprintf(outfile,'\n*** MODEL DESCRIPTION ***\n');
    fprintf(outfile,'.tran 10m\n');
    fclose(outfile);
    s=0;
end
