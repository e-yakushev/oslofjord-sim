FROM europe-west1-docker.pkg.dev/niva-cd/images/fjordsim-oceananigans:c1bc5a4

WORKDIR /app

COPY simulation.jl ./simulation.jl
COPY Oxydep.jl ./Oxydep.jl

ENV SIMULATION_LAUNCHER=/app/simulation.jl
ENV JULIA_DEPOT_PATH=/usr/local/julia_depot
ENV JULIA_PROJECT=/app

COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

ENTRYPOINT ["./entrypoint.sh"]