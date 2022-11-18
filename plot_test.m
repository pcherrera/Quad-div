%% U plot
    
    figure;  
    trisurf(elements,coordinates(:,1),coordinates(:,2),U_1(coordinates));
    title('Exact solution U1(:,1)');
    
    figure;  
    trisurf(elements,coordinates(:,1),coordinates(:,2),U_2(coordinates));
    title('Exact solution U2(:,1)');
    
    figure;
    myP0plot(coordinates,elements,x(1:2:dimU1),pU1);
    title('Approximated solution U1(:,1) ');
    
    figure;
    myP0plot(coordinates,elements,x(2:2:dimU1),pU1);
    title('Approximated solution U2(:,1) ');
    
%% grad divU plot 

    figure;  
    trisurf(elements,coordinates(:,1),coordinates(:,2),gdivU_1(coordinates));
    title('Exact solution gdivU1(:,1)');
    
    figure;  
    trisurf(elements,coordinates(:,1),coordinates(:,2),gdivU_2(coordinates));
    title('Exact solution gdivU2(:,1)');
    
    figure;
    myP0plot(coordinates,elements,x(dimU1+1:2:dimU1+dimU2),pU2);
    title('Approximated solution gdivU1(:,1) ');
    
    figure;
    myP0plot(coordinates,elements,x(dimU1+2:2:dimU1+dimU2),pU2);
    title('Approximated solution gdivU2(:,1) ');
        
    
   %% Trace approximation 
   figure
   trisurf(elements,coordinates(:,1),coordinates(:,2),x(dimU+(1:2:nC(j))'))
   title('Trace approximation');
    
   %% relative error 
   figure 
   loglog(dofDPG,[rel_errU1;rel_errU1P],'LineWidth',2.5)
   title('Relative error %')
   legend('Approximated solution','L2-Projection')